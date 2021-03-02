from fpdf import FPDF
import os


for i in range(2,7):
    pdf = FPDF()
    for image in sorted(os.listdir("R=0.%d"%i)):
        if "tree" in image or "scatter" in image:
            continue
        pdf.add_page()
        pdf.set_font('Arial', 'B', 16)
        pdf.cell(40, 50, image)
        print(image+"R=0.%d"%i)
        pdf.image("R=0.%d/"%i+image, 10,60,200,180)
    pdf.output('R=0.%i.pdf'%i, 'F')