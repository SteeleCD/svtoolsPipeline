import sys
args = sys.argv

# arguments
cn_file = args[1]
fusion_file=args[2]
chrom=args[3]
outFile=args[4]

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
    # The arguments below are optional:
    xlabel = "position (Mb)",
    logbase = 4,
    ymin = 0.5,
    ymax = 200,
    yticks = [1,4,16,64]
    )

