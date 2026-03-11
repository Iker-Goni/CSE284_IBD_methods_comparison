import gzip

vcf_file = "../data/lwk_phased_clean.vcf.gz"
ped_file = "../data/lwk_phased.ped"
map_file = "../data/lwk_phased.map"

samples = []
genotypes = []

with gzip.open(vcf_file, "rt") as f:
    for line in f:
        if line.startswith("##"):
            continue
        
        if line.startswith("#CHROM"):
            parts = line.strip().split()
            samples = parts[9:]
            genotypes = [[] for _ in samples]
            continue

        parts = line.strip().split()

        chrom = parts[0]
        snp = parts[2]
        pos = parts[1]
        ref = parts[3]
        alt = parts[4]

        # write MAP line
        with open(map_file, "a") as m:
            m.write(f"{chrom}\t{snp}\t0\t{pos}\n")

        gts = parts[9:]

        for i, gt in enumerate(gts):
            gt = gt.split(":")[0]

            if "|" in gt:
                a, b = gt.split("|")
            else:
                a, b = gt.split("/")

            a = ref if a == "0" else alt
            b = ref if b == "0" else alt

            genotypes[i].extend([a, b])

with open(ped_file, "w") as p:
    for sample, gts in zip(samples, genotypes):
        row = [sample, sample, "0", "0", "0", "-9"] + gts
        p.write(" ".join(row) + "\n")
