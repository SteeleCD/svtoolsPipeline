import sys
args = sys.argv

# arguments
fusion_file=args[1]
chrom=args[2]
outDir=args[3]

# read fusions
fusions = sv_data.get_fusions(fusion_file, chrom)

# kc_tests
kc_tests.test_E1(
            fusions,
            outDir+"kcE1.png"
            )


