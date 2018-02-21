import sys
args = sys.argv

# arguments
cn_file = args[1]
fusion_file=args[2]
chrom=args[3]
outFile=args[4]
sample=args[5]

if chrom=="chrX":
    vert = [77]
    gene=["ATRX"]
elif chrom=="chr17":
    vert = [7.57]
    gene=["TP53"]
elif chrom=="chr5":
    vert = [1.27]
    gene=["TERT"]
else:
    vert=[-100]
    gene="none"

# import sv tools
from sv_tools import sv_data, sv_diagram


# read CN
x, cn = sv_data.get_x_cn(cn_file, chrom)
# read fusions
fusions = sv_data.get_fusions(fusion_file, chrom)

# sv diagram
sv_diagram.plot_sv_diagram(
    x, cn, fusions,
    outfile = outFile,
    vert=vert,
    genes=gene,
    chrom=chrom,
    sample=sample,
    # The arguments below are optional:
    xlabel = "position (Mb)",
    logbase = 4,
    ymin = 0.5,
    ymax = 200,
    yticks = [1,4,16,64]
    )

