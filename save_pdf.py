import time
import data
from decoot import Output
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, 
                               Image, Table, TableStyle, PageBreak)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
import datetime
import io
import os



 
def save_pdf(output, name):
    styles=getSampleStyleSheet()
    now = datetime.datetime.now()
    doc = SimpleDocTemplate(name + ".pdf",pagesize=letter,
                            rightMargin=72,leftMargin=72,
                            topMargin=72,bottomMargin=18)
    Report=[]
    logo = "decoot.png"
    formatted_time = time.ctime()
    im = Image(logo, 2*inch, 1*inch)
    Report.append(im)
    Report.append(Spacer(1, 12))
    styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))
    ptext = '<font size=12>%s</font>' % formatted_time
    Report.append(Paragraph(ptext, styles["Normal"]))
    Report.append(Spacer(1, 12))



    imgByteArr = io.BytesIO()
    output.img.save(imgByteArr, format='PNG')
    #imgByteArr = imgByteArr.getvalue()
    h, w = output.img.size
    img = Image(imgByteArr, 8*inch, 8*w/h*inch)
    Report.append(img)
    Report.append(PageBreak())

    imgByteArr = io.BytesIO()
    output.graph_error.save(imgByteArr, format='PNG')
    #imgByteArr = imgByteArr.getvalue()
    h, w = output.img.size
    img = Image(imgByteArr, 8*inch, 8*w/h*inch)
    Report.append(img)
    Report.append(PageBreak())
    
    table = [['', 'expected', 'reached', 'difference']]

    for i, j, k in zip(output.vec2fit, output.reached_distribution, range(21)):
        row = []
        row.append(data.names[k])
        if(i):
            row.append("%.3f" %(i))
        else:
            row.append('0')
        if(j):
            row.append("%.3f      " %(j/sum(output.reached_distribution)))
        else:
            row.append('0')
        if(abs(i-j/sum(output.reached_distribution)) > 10**-3):
            if(i-j/sum(output.reached_distribution) < 0):
                row.append("  %.3f" %(round(-(i-j/sum(output.reached_distribution)), 3)))
            else:
                row.append(" %.3f" %(round(-(i-j/sum(output.reached_distribution)), 3)))

        else:
            row.append('  0')
        table.append(row)

    t = Table(table)

    Report.append(t)

    
    ptext = '<font size=12>Length of amino acid sequence: %d</font>' %output.length 
    Report.append(Paragraph(ptext, styles["Justify"]))
    Report.append(Spacer(1, 12))

    if(output.model_distribution_name):
        ptext = '<font size=12>Codonstring was optimized to be close to %s codon distribution</font>' %output.model_distribution
        Report.append(Paragraph(ptext, styles["Justify"]))
        Report.append(Spacer(1, 12))

    ptext = '<font size=12>Maximum ratio of an amino acid in a single codon: %.2f</font>' %output.threshold 
    Report.append(Paragraph(ptext, styles["Justify"]))
    Report.append(Spacer(1, 12))

    if(output.spiked_codons):
        ptext = '<font size=12>Spiked codons</font>' 
    else:
        ptext = '<font size=12>Degenerate codons</font>' 
    Report.append(Paragraph(ptext, styles["Justify"]))
    Report.append(Spacer(1, 12))


    ptext = '<font size=12>Statistics (in sense of sum of squares of differences)</font>' 
    Report.append(Paragraph(ptext, styles["Justify"]))
    Report.append(Spacer(1, 12))
    
    
    
    ptext = '<font size=12>Mean error of a single protein : %.6f</font>' %output.mean
    Report.append(Paragraph(ptext, styles["Justify"]))
    Report.append(Spacer(1, 12))
    
    
    ptext = '<font size=12>Variance of errors of a single protein : %.6f</font>' %output.var
    Report.append(Paragraph(ptext, styles["Justify"]))
    Report.append(Spacer(1, 12))

    ptext = '<font size=12>GC content : %.6f</font>' %output.gc
    Report.append(Paragraph(ptext, styles["Justify"]))
    Report.append(Spacer(1, 12))


    ptext = '<font size=12>weight : %.6f</font>' %output.weight
    Report.append(Paragraph(ptext, styles["Justify"]))
    Report.append(Spacer(1, 12))


    ptext = '<font size=12> Removed triplets: '
    if(len(output.removed_triplets)):
        for triplet in output.removed_triplets:
            ptext += triplet + " "
    else:
        ptext += "None"
    ptext += "</font>" 
    Report.append(Paragraph(ptext, styles["Justify"]))
    Report.append(Spacer(1, 12))
    
    ptext = '<font size=12>Nucleotide sequence in text format:</font>' 
    Report.append(Paragraph(ptext, styles["Justify"]))
    Report.append(Spacer(1, 12))

    ptext = output.output_string
    Report.append(Paragraph(ptext, styles["Justify"]))
    Report.append(Spacer(1, 12))

    doc.build(Report)