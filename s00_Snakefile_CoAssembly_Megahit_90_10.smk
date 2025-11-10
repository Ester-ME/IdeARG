configfile: '/mnt/thor/bigdata/ideARG_2020/results_coassembly_megahit_MAGs_90_10/s00_config_CoAssembly_Megahit_90_10.yaml'


rule all:
    input:
        expand("{wdir}/{sample}/prodigal", sample=config["SAMPLES"], wdir=config["WDIR"]),
        expand("{wdir}/{sample}/ASSEMBLY_hmms_AMRFinderPlus_CONSERVEDMATCHES_cutGA", sample=config["SAMPLES"], wdir=config["WDIR"]),
        expand("{wdir}/{sample}/mobileOG", wdir=config["WDIR"], sample=config["SAMPLES"]),
        expand("{wdir}/{sample}/diamond_heavyMetalResistance", wdir=config["WDIR"], sample=config["SAMPLES"]),
        expand("{wdir}/{sample}/mags_prodigal", sample=config["SAMPLES"], wdir=config["WDIR"]),
        expand("{wdir}/{sample}/mags_hmms_AMRFinderPlus_CONSERVEDMATCHES_cutGA", sample=config["SAMPLES"], wdir=config["WDIR"]),
        expand("{wdir}/{sample}/mags_mobileOG", wdir=config["WDIR"], sample=config["SAMPLES"]),



#####
#####
##### ASSEMBLY - BASED ANALYSIS
#####
#####

# ## GeoMosaic conda environment
rule run_prodigal:
    input:
        contig_path="{wdir}/{sample}/megahit",
    output:
        directory("{wdir}/{sample}/prodigal")
    params:
        extra="-p meta",
        quiet="-q"
    threads: 1
    run:
        shell("mkdir -p {output}")
        shell("prodigal -i {input.contig_path}/contigs.fasta -o {output}/genes.out -a {output}/protein_translations.faa  {params.quiet} {params.extra}")

        from geomosaic.parser.prodigal_orf_mapping import parsing_prodigal_orfs

        fasta_input = f"{output}/protein_translations.faa"
        output_mapping = f"{output}/orf_contig_mapping.tsv"
        output_fasta = f"{output}/orf_predicted.faa"
        output_simple_mapping = f"{output}/simple_orf_contig_mapping.tsv"

        parsing_prodigal_orfs(fasta_input, output_mapping, output_fasta, output_simple_mapping)


# conda activate mobileOG
rule run_DIAMOND_HeavyMetalResistance:
    input:
        prodigal_path="{wdir}/{sample}/prodigal",
        hmr_db="/mnt/thor/bigdata/ideARG_2020/heavy_metal_resistance/args_oap_metaldb/heavyMetalResistance.dmnd"
    output:
        folder=directory("{wdir}/{sample}/diamond_heavyMetalResistance")
    threads: 5
    run:    
        shell("mkdir -p {output.folder}")
        
        shell("diamond blastp -d {input.hmr_db} \
                -q {input.prodigal_path}/orf_predicted.faa \
                -o {output.folder}/{wildcards.sample}.hmr.tsv \
                --sensitive \
                --id 60 \
                -e 1e-10 \
                --threads {threads} \
                -k 20 \
                --outfmt 6 sseqid qseqid stitle qtitle pident bitscore evalue length slen qlen sstart send qstart qend qcovhsp scovhsp mismatch gapopen")



# conda activate mobileOG
rule run_mobileOG:
    input:
        prodigal_path="{wdir}/{sample}/prodigal",
        mobileOG_db="/mnt/thor/bigdata/ideARG_2020/mobileOG_database/mobileOG-db_beatrix-1.6.dmnd"
    output:
        folder=directory("{wdir}/{sample}/mobileOG")
    threads: 5
    run:    
        shell("mkdir -p {output.folder}")
        
        shell("diamond blastp -d {input.mobileOG_db} \
                -q {input.prodigal_path}/orf_predicted.faa \
                -o {output.folder}/{wildcards.sample}.mobileog.m8 \
                --sensitive \
                --id 60 \
                -e 1e-10 \
                -k 20 \
                --outfmt 6 sseqid qseqid stitle qtitle pident bitscore evalue length slen qlen sstart send qstart qend qcovhsp scovhsp mismatch gapopen")


rule run_ASSEMBLY_hmmsearch_AMRFPlus_conserved_matches_cutGA:
    input:
        prodigal_path="{wdir}/{sample}/prodigal",
    output:
        folder=directory("{wdir}/{sample}/ASSEMBLY_hmms_AMRFinderPlus_CONSERVEDMATCHES_cutGA")
    params:
        hmm_folder=config["hmm_folder_AMRFinderPlus"],
        local_sample="{sample}"
    threads: 
        5
    run:
        shell("mkdir -p {output.folder}")
        import pandas as pd

        df_mapping = pd.read_csv(os.path.join(str(input.prodigal_path), "simple_orf_contig_mapping.tsv"), sep="\t")

        results_filename = os.path.join(str(output.folder), "presence_results.tsv")

        list_output_files = []
        for hmm in os.listdir(params.hmm_folder):
            if not hmm.endswith('.HMM'):
                continue
            
            filename = hmm.split(".HMM")[0]
            out_path = os.path.join(str(output.folder), "all_models_results", filename)
            shell("mkdir -p {out_path}")

            hmm_file = os.path.join(params.hmm_folder, hmm)
            shell("hmmsearch --tblout /dev/null -o {out_path}/hmmsearch_output.txt --cut_ga --cpu {threads} --notextw {hmm_file} {input.prodigal_path}/orf_predicted.faa") 
            if os.stat(os.path.join(out_path, "hmmsearch_output.txt")).st_size == 0:
                continue

            list_output_files.append(str(os.path.join(out_path, "hmmsearch_output.txt")))
            
        from geomosaic.parser.make_hmmsearch_dataframe import make_hmmsearch_dataframe
        df_hmmresults = make_hmmsearch_dataframe(list_output_files, mags=False)
        df_hmmresults.drop_duplicates(inplace=True)

        m1 = df_hmmresults.merge(df_mapping, on="orf_id", how="left")
        m1["sample"] = str(params.local_sample)

        m1.to_csv(results_filename, sep="\t", header=True, index=False)

        shell("( cd {output.folder} && rm -r all_models_results )")



#####
#####
##### MAGS - BASED ANALYSIS
#####
#####

rule run_mags_prodigal:
    input:
        mags_folder=expand("{wdir}/{sample}/MAGs", allow_missing=True),
    output:
        directory("{wdir}/{sample}/mags_prodigal")
    params:
        extra="-p meta",
        quiet="-q"
    threads: 1
    run:
        shell("mkdir -p {output}")
        
        import pandas as pd
        from geomosaic.parser.prodigal_orf_mapping import parsing_prodigal_orfs_MAGs

        df_mags = pd.read_csv(str(os.path.join(str(input.mags_folder), "MAGs.tsv")), sep="\t")
        
        mags_list_file = os.path.join(str(output), "mags_list.txt")
        with open(mags_list_file, "wt") as fd:
            for mag in list(df_mags.MAGs):
                fd.write(f"{mag}\n")
                print(f"Computing ORF prediction for {mag}")
                output_folder_mag=str(os.path.join(str(output), mag))

                shell("mkdir -p {output_folder_mag}")
                shell("prodigal -i {input.mags_folder}/fasta/{mag}.fasta \
                        -o {output_folder_mag}/genes.gff \
                        -a {output_folder_mag}/protein_translations.faa \
                        -f gff \
                        {params.extra} \
                        {params.quiet}")

                fasta_input = str(os.path.join(output_folder_mag, "protein_translations.faa"))
                output_mapping = str(os.path.join(output_folder_mag, "orf_contig_mapping.tsv"))
                output_fasta = str(os.path.join(output_folder_mag, "orf_predicted.faa"))
                output_simple_mapping = str(os.path.join(output_folder_mag, "simple_orf_contig_mapping.tsv"))

                parsing_prodigal_orfs_MAGs(fasta_input, output_mapping, output_fasta, output_simple_mapping)


rule run_mags_hmmsearch_AMRFPlus_conserved_matches_cutGA:
    input:
        prodigal_path="{wdir}/{sample}/mags_prodigal",
    output:
        folder=directory("{wdir}/{sample}/mags_hmms_AMRFinderPlus_CONSERVEDMATCHES_cutGA")
    params:
        hmm_folder=config["hmm_folder_AMRFinderPlus"],
        local_sample="{sample}"
    threads: 
        5
    run:
        shell("mkdir -p {output.folder}")
        import pandas as pd

        mags_list = []
        with open(os.path.join(str(input.prodigal_path), "mags_list.txt")) as magslist_fd:
            for x in magslist_fd:
                mags_list.append(x.rstrip("\n"))

        for mag in mags_list:
            mag_orf_folder = os.path.join(str(input.prodigal_path), mag)

            df_mapping = pd.read_csv(os.path.join(mag_orf_folder, "simple_orf_contig_mapping.tsv"), sep="\t")

            mags_temp_outfolder = os.path.join(str(output.folder), mag)
            results_filename = os.path.join(str(output.folder), f"{mag}_presence_results.tsv")

            list_output_files = []
            for hmm in os.listdir(params.hmm_folder):
                if not hmm.endswith('.HMM'):
                    continue
                
                filename = hmm.split(".HMM")[0]
                out_path = os.path.join(mags_temp_outfolder, filename)
                shell("mkdir -p {out_path}")

                hmm_file = os.path.join(params.hmm_folder, hmm)
                shell("hmmsearch --tblout /dev/null -o {out_path}/hmmsearch_output.txt --cut_ga --cpu {threads} --notextw {hmm_file} {mag_orf_folder}/orf_predicted.faa") 
                if os.stat(os.path.join(out_path, "hmmsearch_output.txt")).st_size == 0:
                    continue

                list_output_files.append(str(os.path.join(out_path, "hmmsearch_output.txt")))
            
            from geomosaic.parser.make_hmmsearch_dataframe import make_hmmsearch_dataframe
            df_hmmresults = make_hmmsearch_dataframe(list_output_files, mags=True)
            df_hmmresults.drop_duplicates(inplace=True)

            m1 = df_hmmresults.merge(df_mapping, on="orf_id", how="left")
            m1["sample"] = str(params.local_sample)

            m1.to_csv(results_filename, sep="\t", header=True, index=False)

            shell("(cd {output.folder} && rm -r {mag})")


# conda activate mobileOG
rule run_mags_mobileOG:
    input:
        prodigal_path="{wdir}/{sample}/mags_prodigal",
        mobileOG_db="/mnt/thor/bigdata/ideARG_2020/mobileOG_database/mobileOG-db_beatrix-1.6.dmnd"
    output:
        folder=directory("{wdir}/{sample}/mags_mobileOG")
    threads: 5
    run:
        mags_list = []
        with open(os.path.join(str(input.prodigal_path), "mags_list.txt")) as magslist_fd:
            for x in magslist_fd:
                mags_list.append(x.rstrip("\n"))

        for mag in mags_list:
            mag_orf_predicted = os.path.join(str(input.prodigal_path), mag, "orf_predicted.faa")
            mags_outfolder = os.path.join(str(output.folder), mag)
            
            shell("mkdir -p {mags_outfolder}")
            
            shell("diamond blastp -d {input.mobileOG_db} \
                    -q {mag_orf_predicted} \
                    -o {mags_outfolder}/{mag}.mobileog.m8 \
                    --sensitive \
                    --id 60 \
                    -e 1e-10 \
                    -k 20 \
                    --outfmt 6 sseqid qseqid stitle qtitle pident bitscore evalue length slen qlen sstart send qstart qend qcovhsp scovhsp mismatch gapopen")

