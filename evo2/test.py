import pysam

vcf_path = "./benchmark/ClinVar_Coding_SNV_PB.vcf"
with pysam.VariantFile(vcf_path) as vcf:
    # Print header INFO fields
    print("Available INFO fields:")
    for field in vcf.header.info:
        print(f"- {field}: {vcf.header.info[field]}")
    
    # Print first 5 records with their INFO fields
    print("\nFirst 5 records:")
    for i, record in enumerate(vcf):
        if i >= 5:
            break
        print(f"\nRecord {i+1}:")
        print(f"CHROM: {record.chrom}, POS: {record.pos}")
        print("INFO fields:")
        for field in record.header.info:
            if field in record.info:
                value = record.info[field]
                print(f"- {field}: {value} (type: {type(value)})")