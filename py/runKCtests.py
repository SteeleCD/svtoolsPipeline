import sys
args = sys.argv

# arguments
fusion_file=args[1]
chrom=args[2]
outDir=args[3]

from sv_tools import sv_data, kc_tests

# read fusions
fusions = sv_data.get_fusions(fusion_file, chrom)

# kc_tests
#Test A - Clustering of breakpoints (exponential QQ plot)
#Test B - Oscillation among a few CN states
#Test C - Interspersed loss and retention of heterozygosity
#Test D - Rearrangements affecting only one haplotype
#Test E1 - Randomness of fragment joins
#Test E2 - Randomness of fragment order
#Test F - Ability to walk the derivative chromosome.
#A-D & E2 not implemented in sv_tools
# E1 
kc_tests.test_E1(
            fusions,
            outDir+"kcE1.png"
            )
# F - not sure how to define walk from data (made of H/T combinations)
#kc_tests.F_walk(
#            "HTHTTHTHTH"
#            )

