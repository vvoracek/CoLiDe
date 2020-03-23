import data
import math
from copy import deepcopy as dc
from random import choice, random, shuffle as randshuffle
from itertools import product
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import io
from Bio.SeqUtils.CodonUsage import SynonymousCodons
from Bio.Data.IUPACData import protein_letters_1to3
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import translate
from Bio.SeqUtils import GC
from PIL import Image, ImageTk
from sys import stderr



class BasesToCodons:
    """
    Wrapper class for a list of lists for codon translation. 

    For example [['TTT', 'TTC'], ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'], ...]
    """
    def __init__(self):
        stop = ['TAA', 'TAG', 'TGA']
        self.codons = dc([SynonymousCodons[protein_letters_1to3[i.replace(" ", "")].upper()] for i in data.names[:-1]] + [stop])

class Tribase:
    """
    Wrapper class for list of lists containing bases.

    For example [[0,0,100,0], [25, 40, 5, 30], [50,50,0,0]]
    stands for three bases, where the second is 25% T, 40% C, 5% A and 30%G.

    """
    def __init__(self, inarg, b2c):
        if(type(inarg) == list):
            self.bases = inarg
            self.set_tribase_codons(b2c)
        else:
            self.codons = dc(inarg.codons)
            self.bases = dc(inarg.bases)

    def set_tribase_codons(self, b2c):
        transposition_table = {'T' : 0, 
                               'C' : 1, 
                               'A' : 2, 
                               'G' : 3}
        def prob(string):
            ret = self.bases[0][transposition_table[string[0]]]
            ret *= self.bases[1][transposition_table[string[1]]]
            ret *= self.bases[2][transposition_table[string[2]]]
            return ret

        self.codons = [sum(map(prob, i)) for i in b2c.codons]

    def generate_random_tribase(parameters, d = 0):
        if(d == 100):
            return None
        bases = dc(choice(parameters.valid_tribases))
        
        for i in range(3):
            tmp = []
            while(1):
                if(parameters.spiked_codons):
                    tmp = [random() for _ in range(sum(bases[i]))]
                    tmp = [round(100*i/sum(tmp)) for i in tmp]
                else:
                    tmp = [100//sum(bases[i]) for _ in range(sum(bases[i]))]
                    if(sum(tmp) != 100):
                        tmp[int(random()*3)] += 1
                    
                if(sum(tmp) == 100):
                    break
            ptr1 = ptr2 = 0
            while(ptr1<4):
                if(bases[i][ptr1]):
                    bases[i][ptr1] = tmp[ptr2]
                    ptr2 += 1
                ptr1 += 1
            
        ret = Tribase(bases,parameters.b2c)

        if(max(ret.codons) > parameters.threshold * 10**6):
            return Tribase.generate_random_tribase(parameters, d+1)
        return ret
    
    def generate_subtribases(self):
        ret = []
        idxs = [[idx for idx, element in enumerate(base) if element != 0] for base in self.bases]
        for i, j, k in product(*idxs):
            to_add = [[0,0,0,0], [0,0,0,0], [0,0,0,0]]
            to_add[0][i] = 1
            to_add[1][j] = 1
            to_add[2][k] = 1
            ret.append(to_add)
        return ret


class TribaseString:
    def __init__(self, parameters):
        self.parameters = parameters
        self.tribases = parameters.length * [0]
        self.codons = None
        self.combined_bases = 64*[0]
        self.failed_init = 0

        for i in range(parameters.length):
            self.tribases[i] = Tribase.generate_random_tribase(parameters)
            if(self.tribases[i] == None):
                self.failed_init = 1
                return
    
        self.set_codons()
        self.set_combined_bases()

        self.sum_of_codons = sum(self.codons)
        self.compute_error()
        if(self.parameters.model_distribution):
            self.sum_of_bases = sum(self.combined_bases)
            self.compute_error_model_distribution()

    def set_codons(self):
        self.codons = len(self.tribases[0].codons) * [0]
        for i in range(self.parameters.length):
            for j in range(len(self.tribases[i].codons)):
                self.codons[j] += self.tribases[i].codons[j]
                
    def set_combined_bases(self):
        self.combined_bases = 64*[0]
        for i in range(self.parameters.length):
            for a in range(4):
                for b in range(4):
                    for c in range(4):
                        self.combined_bases[16*a + 4*b + c] += self.tribases[i].bases[0][a] *\
                                                            self.tribases[i].bases[1][b] *\
                                                            self.tribases[i].bases[2][c]

    def compute_error(self):
        self.error = sum([(self.codons[i]/self.sum_of_codons-self.parameters.vec2fit[i])**2 for i in range(len(self.parameters.vec2fit))])


    def compute_error_model_distribution(self):
        self.error_model_distribution =  sum([(self.combined_bases[i]/self.sum_of_bases-self.parameters.model_distribution[i])**2 for i in range(len(self.parameters.model_distribution))])


    def update(self, optimize_model_distribution = False):
        random_index = int(random()*self.parameters.length)
        new_tribase = Tribase.generate_random_tribase(self.parameters)
        if(not new_tribase or max(new_tribase.codons) >  self.parameters.threshold * 10**6):
            return False

        if(optimize_model_distribution):
            bas  = dc(self.combined_bases)
            sum_of_bases = self.sum_of_bases

            for i in range(len(bas)):
                bas[i] = bas[i] + new_tribase.bases[0][i//16] *\
                                  new_tribase.bases[1][(i%16)//4] *\
                                  new_tribase.bases[2][i%4]\
                                - self.tribases[random_index].bases[0][i//16] *\
                                  self.tribases[random_index].bases[1][(i%16)//4] *\
                                  self.tribases[random_index].bases[2][i%4]
            new_error = sum([(bas[i]/sum_of_bases-self.parameters.model_distribution[i])**2\
                 for i in range(len(self.parameters.model_distribution))])

            if(self.error_model_distribution > new_error):
                self.error_model_distribution = new_error
                self.combined_bases = bas
                self.tribases[random_index] = new_tribase
                return True
        else:
            cod = dc(self.codons)
            cod1 = new_tribase.codons
            cod2 = self.tribases[random_index].codons
            
            for i in range(len(cod)):
                cod[i] = cod[i] + cod1[i] - cod2[i]
            new_error = sum([(cod[i]/self.sum_of_codons-self.parameters.vec2fit[i])**2 for i in range(len(self.parameters.vec2fit))])

            if(self.error > new_error):
                self.error = new_error
                self.codons = cod
                self.tribases[random_index] = new_tribase
                return True
        return False

    def sample_protein(self):
        codons = len(self.parameters.b2c.codons) *[0]
        code = ""
        for tribase in self.tribases:
            bases = tribase.bases
            codon = [[0,0,0,0], [0,0,0,0], [0,0,0,0]]
            for i in range(len(bases)):
                base = bases[i]
                r = int(100*random())+1
                cumsum = 0
                for j in range(len(base)):
                    cumsum += base[j]
                    if(cumsum >= r):
                        codon[i][j] = 1
                        break

            t = Tribase(codon, self.parameters.b2c)
            code += translate_triplets(codon)
            codons = [i + j for i, j in zip(codons, t.codons)]
        PA = ProteinAnalysis(translate(code))
        gc = GC(code)
        try:
            w = PA.molecular_weight()
        except:
            w = 0
        return codons, gc, w

    def get_statistics(self, number_of_samples = 10**3):
        errors = []
        weights = []
        gcs = []
        for i in range(number_of_samples):
            codons, gc, weight = self.sample_protein()
            gcs.append(gc)
            weights.append(weight)
            error = (sum([(codons[i]/sum(codons)-self.parameters.vec2fit[i])**2 for i in range(len(self.parameters.vec2fit))]))
            errors.append(error)

        mean = sum(errors)/len(errors)
        var = sum((error - mean)**2 for error in errors)/(len(errors)-1)

        return mean, var, sum(gcs)/len(gcs), sum(weights)/len(weights)


def translate_triplets(argin):
    transposition_table = {'T' : 0, 
                            0  :'T',
                           'C' : 1,
                            1  :'C',
                           'A' : 2, 
                            2  :'A',
                           'G' : 3,
                            3  :'G'}
    if(type(argin) == str):
        ret = [[0,0,0,0], [0,0,0,0], [0,0,0,0]]
        for i, j in zip(string, range(len(string))):
            ret[j][transposition_table[i]] = 1
        return ret
    else:
        ret = ""
        for i in range(3):
            for j in range(4):
                if(argin[i][j]):
                    ret += transposition_table[j]
        if(len(ret) == 3):
            return ret
        else:
            print("translate_triplets", ret)




def subset(b1, b2):
    for i,j in zip(b1,b2):
        for k, l in zip(i,j):
            if(k > l):
                return False
    return True


class Parameters:
    def __init__(self, model_distribution, threshold, spiked_codons, removed_triplets, vec2fit, length, base2codon):
        self.model_distribution = model_distribution
        self.threshold = threshold
        self.spiked_codons = spiked_codons
        self.removed_triplets = removed_triplets
        self.vec2fit = vec2fit
        self.length = length
        self.b2c = base2codon
    
        self.removed_codons = [j for i, j in zip(vec2fit, range(len(vec2fit))) if(not i)]
        self.valid_tribases = []

    def __repr__(self):
        return "model distribution: " + str(self.model_distribution) +"\n"+\
        "threshold: " + str(self.threshold)+ "\n"+\
        "spiked_codons: " + str(self.spiked_codons)+"\n"+\
        "removed_triplets: " + str(self.removed_triplets)+"\n"+\
        "vec2fit: " + str(self.vec2fit)+"\n"+\
        "length: " + str(self.length)

class Output:
    def __init__(self, parameters, tribase_string):
        self.tribase_string = tribase_string
        self.parameters = parameters
        self.log_entropy = 0
        self.valid = 1
        for tribase in tribase_string.tribases:
            self.log_entropy += math.log(sum([i*math.log(i) for i in tribase.codons if (i)]))

        self.parameters = parameters
        self.removed_triplets = parameters.removed_triplets
        self.vec2fit = parameters.vec2fit
        self.reached_distribution = tribase_string.codons
        self.mean, self.var, self.gc, self.weight= tribase_string.get_statistics()
        self.spiked_codons = parameters.spiked_codons
        self.length = parameters.length
        self.threshold = parameters.threshold
        self.output_string = Output.get_output_string(
                        tribase_string.tribases, parameters.spiked_codons)
        self.model_distribution_name = Output.get_model_distribution_name(
                        parameters.model_distribution)
        self.graph_error = None
        self.img = None

    def make_imgs(self):
        self.graph_error = self.make_graph_error()

        self.img = Output.make_img(self.tribase_string.tribases)

    def shuffle(self):
        randshuffle(self.tribase_string.tribases)
        self.img = Output.make_img(self.tribase_string.tribases)
        self.output_string = Output.get_output_string(
            self.tribase_string.tribases, self.spiked_codons)

    def get_model_distribution_name(model_distribution):
        model_distribution_name = ""
        for name in data.options:
            if(model_distribution == data.options[name]):
                model_distribution_name = name
        return model_distribution_name

    def get_output_string(tribases, spiked):
        d = {}
        d[1]  = 'T'
        d[2]  = 'C'
        d[3]  = 'Y' # Y = TC
        d[4]  = 'A' 
        d[5]  = 'W' # W = AT
        d[6]  = 'M' # M = AC
        d[7]  = 'H' # H = ACT 
        d[8]  = 'G'
        d[9]  = 'K' # K = GT
        d[10] = 'S' # S = GC
        d[11] = 'B' # B = GCT
        d[12] = 'R' # R = GA
        d[13] = 'D' # D = GAT
        d[14] = 'V' # V = GAC
        d[15] = 'N' # N = GACT
        
        s = ""
        if(not spiked):
            for i in tribases:
                bases = i.bases
                for position in bases:
                    num = 0
                    for k in position[::-1]:
                        num *= 2
                        num += int(bool(k))
                    s += d[num]
            
        else:
            tmp = 16*[1]
            for i in tribases:
                i = i.bases
                for j in i:
                    num = 0
                    for k in j[::-1]:
                        num *= 2
                        num += int(bool(k))
                    if(sum([int(i) for i in bin(num)[2:]]) == 1):
                        s += str(d[num])
                    else:
                        s += "("+str(d[num])+str(tmp[num])+":"
                        s += str(100+round(j[2]))[1:]
                        s += str(100+round(j[1]))[1:]
                        s += str(100+round(j[3]))[1:]
                        s += str(100+round(j[0]))[1:]
                        s += ")"
                        tmp[num] += 1
        return s

    def make_img(tribases):
        imgs =  []
        for idx, codon in zip(range(len(tribases)), tribases):
            sizes = codon.codons
            while(len(sizes) > len(data.colors)):
                data.colors.append((1,0,1))
            colors = [data.colors[i] for i in range(len(sizes)) if(sizes[i])]

            plt.figure(figsize=(4,4))
            plt.pie(sizes, colors=data.colors, shadow=True)
            plt.title(str(idx+1), fontsize = 30)

            buf = io.BytesIO()
            plt.savefig(buf, format='png')
            buf.seek(0)
            imgs.append(buf)
            plt.clf()
            plt.close()
        return Output.concat_images(imgs)

    def concat_images(imgs):
        images = list(map(Image.open, imgs))
        a = int(len(imgs)**0.5)
        b = int(len(imgs)/a+0.99)

        w, h = images[0].size
        
        result = Image.new("RGBA", (b*w, a*h), color='white')
        x = 0
        idx = 0
        for i in range(a):
            for j in range(b):
                if(idx >= len(images)):
                    break
                result.paste(images[idx], (j*w, i*h))
                idx += 1
        return result
    
    def make_graph_error(self):
        N = len(self.parameters.b2c.codons)
        ind = np.arange(N)
        width = 0.3      

        fig = plt.figure()
        ax = fig.add_subplot(111)


        expected = self.vec2fit
        rects1 = ax.bar(ind+width, expected, width, color='black')
        reached = list(map(lambda x:x/sum(self.reached_distribution), self.reached_distribution))
        rects2 = ax.bar(ind, reached, width, color='r')

        ax.set_xticks(ind+width/2)
        ax.set_xticklabels(data.names[:-1] + [" * "] + sum(self.parameters.b2c.codons[21:], []), rotation='vertical')
        #ax.set_xticklabels(dc(data.names), rotation='vertical')

        ax.legend( (rects1[0], rects2[0]), ('expected', 'reached') )

        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        return Image.open(buf)

class Decoot:
    def __init__(self, parameters):
        self.parameters = parameters

    def set_model_distribution(self):
        for val, idx in zip(self.parameters.vec2fit, range(len(self.parameters.vec2fit))):
            for i in range(4):
                for j in range(4):
                    for k in range(4): 
                        l = [[0,0,0,0,], [0,0,0,0], [0,0,0,0]]
                        l[0][i] = 1
                        l[1][j] = 1
                        l[2][k] = 1
                        
                        t = Tribase(l, self.parameters.b2c)
                        
                        if(t.codons[idx]):
                            self.parameters.model_distribution[16*i+4*j+k] *= val
        self.parameters.model_distribution = [i/sum(self.parameters.model_distribution) for i in self.parameters.model_distribution]

    def generate_valid_tribases(self):
        for tmp in range((2**4)**3):
            tmp += (2**4)**3
            tmp = bin(tmp)[3:]
            to_add = [[int(i) for i in tmp[4*it:4*it+4]]for it in range(3)]
            tribase = Tribase(to_add, self.parameters.b2c)
            
            if(sum([bool(i) for i in tribase.codons]) < 1 or \
                (sum([int(bool(i)) for i in tribase.codons]) == 1 and self.parameters.threshold != 1)):
                continue
            add = True

            for i in self.parameters.removed_codons:
                if(tribase.codons[i]):
                    add = False
                    break
            for i in self.parameters.removed_triplets:
                if(add and subset(i, to_add)):
                    add = False
                    break
            if(add):
                self.parameters.valid_tribases.append(to_add)
    
    def __call__(self, out = None):
        self.parameters.valid_tribases = []
        self.generate_valid_tribases()
        if(not self.parameters.valid_tribases):
            return None
        if(self.parameters.model_distribution):
            self.set_model_distribution()

        tribase_string = TribaseString(self.parameters)
        if(tribase_string.failed_init):
            return None


        not_updated = 0
        if(self.parameters.model_distribution):
            while(not_updated < 10*tribase_string.parameters.length):
                not_updated += 1
                if(tribase_string.update(True)):
                    not_updated = 0


        tribase_string.set_codons()
        tribase_string.compute_error()

        not_updated = 0
        while(not_updated < 1000*tribase_string.parameters.length):
            not_updated += 1
            if(tribase_string.update()):
                not_updated = 0

        out = Output(self.parameters, tribase_string)

        if(out.valid):
            return out
        else:
            return None