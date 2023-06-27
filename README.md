# Salmon Proteogenomics Project #
<br />

Software Requirements
---------------
<details><summary><b>WSL (Ubuntu-20.04)</b></summary>
<p>
<dl><dd>

<!-- WSL -->
```sh
wsl --install -d Ubuntu-20.04
```

<!-- libgomp1 -->
* **libgomp1**
<dl><dd><dl><dd>

```sh
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install libgomp1
```

</dd></dl></dd></dl> 

<!-- Java -->
* **Java**
<dl><dd><dl><dd>

```sh
sudo apt-get update
sudo apt-get install default-jre
```

</dd></dl></dd></dl>  

<!-- Anaconda 3 -->
* **Anaconda 3**
<dl><dd><dl><dd>

```sh
wget https://repo.anaconda.com/archive/Anaconda3-2023.03-1-Linux-x86_64.sh
bash Anaconda3-2023.03-1-Linux-x86_64.sh
```

</dd></dl></dd></dl> 

<!-- biopython -->
* **biopython (Python package)**
<dl><dd><dl><dd>

```sh
pip install biopython
```

</dd></dl></dd></dl> 

<!-- augustus-3.3.3 -->
* **augustus-3.3.3**
<dl><dd><dl><dd>

```sh
wget https://github.com/Gaius-Augustus/Augustus/releases/download/v3.3.3/augustus-3.3.3.tar.gz
tar -xzvf augustus-3.3.3.tar.gz
```

</dd></dl></dd></dl> 

<!-- ncbi-blast-2.14.0+ -->
* **ncbi-blast-2.14.0+**
<dl><dd><dl><dd>

**a. Download ncbi-blast-2.14.0+**
```sh
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.0/ncbi-blast-2.14.0+-x64-linux.tar.gz
tar -xzvf ncbi-blast-2.14.0+-x64-linux.tar.gz
rm ncbi-blast-2.14.0+-x64-linux.tar.gz
```
**b. Add blast database directory to path (also include in .bash_profile)**
```
export BLASTDB=$HOME/Atlantic_Salmon_Proteogenomics_Project/data/reference/Actinopterygii/Actinopterygii_refseq/blast_db
```
**c. Download RefSeq fasta**
ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_other/*.protein.faa.gz 

**d. Prepare blast database**
```
python -u ./src/Prep_blastdb_Actinopterygii_RefSeq.py --input [vertebrate_other] --ouput ./data/reference/Actinopterygii/Actinopterygii_refseq --taxon_list ./data/reference/Actinopterygii/ancestor7898_proteomes.tsv --makeblastdb ./src/ncbi-blast-2.14.0+/bin/makeblastdb
```

</dd></dl></dd></dl>

</dd></dl>
</p>
</details> 


Files
---------------
**1. Download "Atlantic_Salmon_Proteogenomics_Project/data" directory**
[Atlantic_Salmon_Proteogenomics_Project.tar.gz](https://drive.google.com/file/d/1xjPzQTb1HtrY-Zfui1ylwDvVCGLmKgVr/view?usp=sharing)<br />


Pipeline
---------------
>**Note** 
>To enable local Blastx and BlastP runs, this pipeline will use a smaller database representing RefSeq sequences from 106 taxon ids with Actinopterygii (7898) as the ancestor.

<!-- 1. SpliceDB construction -->
### 1. SpliceDB construction ###
<dl><dd>
<!-- option 1 -->
<details><summary><b>[Option 1] <a href="https://proteomics2.ucsd.edu/ProteoSAFe/index.jsp">ProteoSAFe (beta) workflow</a></b></summary>
<p>
<dl><dd>

**a. Salmon Proteogenomics Project | SpliceDB Creation**
![image](https://drive.google.com/uc?export=view&id=1cuaFYhyoAOtGPAXy9knTSe_OqDHCcOyh)
**b. Download SpliceDB Fasta file**
[(Example: heart)](https://proteomics2.ucsd.edu/ProteoSAFe/status.jsp?task=49e26c2650834dad803e15b2973438b3)<br />
![image](https://drive.google.com/uc?export=view&id=1sRDUhteGIDh3A1TZjBNlqanPqsNAuWSS)

</dd></dl>
</p>
</details>

<!-- option 2 -->
<details><summary><b>[Option 2] Local script</b></summary>
<p>
<dl><dd>

**a. Setup conda environment for python 2.7**
```
conda create --name py2 python=2.7
```
**b. Construct SpliceDB**
```
conda activate py2
python -u ./src/SpliceDBcreation/buildSpliceGraph_addSAMonly.py --RAMmemory 10 --dataset 2019_tissues_n23 --min_reads_per_coordinate 3 --dna_fasta_dir ./data/reference/RefSeq/GCF_000233375.1_ICSASG_v2_genomic.fna
conda deactivate
```

</dd></dl>
</p>
</details>
</dd></dl>

<!-- 2. MSGF+ (Known->SpliceDB) -->
### 2. MSGF+ (Known->SpliceDB) ###
<dl><dd>
<!-- a. -->
<details><summary><b>a. Salmon Proteogenomics Project | Comprehensive Database Search (known, splice) <a href="https://proteomics2.ucsd.edu/ProteoSAFe/index.jsp">workflow</a></b></summary>
<p>
<dl><dd>

![image](https://drive.google.com/uc?export=view&id=17gTNJ_Vl_NwSb8Aj7jWHqMkvNt5U-syp)

</dd></dl>
</p>
</details>

<!-- b. -->
<details><summary><b>b. Download output (<a href="https://proteomics2.ucsd.edu/ProteoSAFe/status.jsp?task=1945c8c55b504b7c91cd2cbf9dc5f66b">Example: heart</a>)</b></summary>
<p>
<dl><dd>

![image](https://drive.google.com/uc?export=view&id=1G1kl4TNBmclKvHZVraEN8-pyglf5Yvaf)

</dd></dl>
</p>
</details>

<!-- c. -->
<details><summary><b>c. Placement of "CUSTOM_COMPREHENSIVE_DB_SEARCH-xxxxxxxx-all_events-main.tsv" file</b></summary>
<p>
<dl><dd>

>**Output:**<br />
>./data/**[dataset]**/1_RefSeq_SpliceDB_Search/workflow_output/**[subdir]**/**[downloaded_events_file]**
* **[dataset]** = 2019_tissues_n23<br />
* **[subdir]** = heart<br />
* **[downloaded_events_file]** = CUSTOM_COMPREHENSIVE_DB_SEARCH-xxxxxxxx-all_events-main.tsv<br />

</dd></dl>
</p>
</details>
</dd></dl>

<!-- 3. Construct RefSeq+Ensembl database and run MSGF+ -->
### 3. Construct RefSeq+Ensembl database and run MSGF+ ###
<dl><dd>
<!-- a. -->
<details><summary><b>a. Parse proteogenomic events & construct RefSeq+Ensembl database</b></summary>
<p>
<dl><dd>

```
python -u ./src/Combine_enosi_files.py --dataset 2019_tissues_n23
```

>**Output:**<br />
>./data/**[dataset]**/1_RefSeq_SpliceDB_Search/combined_Enosi_Output/**GCF_000233375.1_ICSASG_v2_protein+Ensembl.fasta**<br />
* **[dataset]** = 2019_tissues_n23<br />

</dd></dl>
</p>
</details>

<!-- b. -->
<details><summary><b>b. Salmon Proteogenomics Project | MSGFPlus_PeptideMatchCMD <a href="https://proteomics2.ucsd.edu/ProteoSAFe/index.jsp">workflow</a></b></summary>
<p>
<dl><dd>

![image](https://drive.google.com/uc?export=view&id=10Cck2wGT6SSD24GkeDlX5yraHKWboKsb)

</dd></dl>
</p>
</details>

<!-- c. -->
<details><summary><b>c. Download output (<a href="https://proteomics2.ucsd.edu/ProteoSAFe/status.jsp?task=13a7eb86c6eb4292bf1cdc66789b17d7">Example: 2019_tissues_n23 dataset</a>)</b></summary>
<p>
<dl><dd>

![image](https://drive.google.com/uc?export=view&id=15EYWm3OtBzjLS-06rkNtuMnwI1H-B80r)

</dd></dl>
</p>
</details>

<!-- d. -->
<details><summary><b>d. Placement of "ExactMatch_output" & "MSGF_combinedtsv" directories</b></summary>
<p>
<dl><dd>

>**Output:**<br />
>./data/**[dataset]**/2_RefSeq_Ensembl_Search/workflow_output/**ExactMatch_output** <br />
>./data/**[dataset]**/2_RefSeq_Ensembl_Search/workflow_output/**MSGF_combinedtsv** <br />

* **[dataset]** = 2019_tissues_n23 

</dd></dl>
</p>
</details>

<!-- e. -->
<details><summary><b>e. Compute FDR and generate "MSGF_fdrdir" directory</b></summary>
<p>
<dl><dd>

```
conda activate py2
python -u ./src/Run_compute_FDR_customDB.py --dataset 2019_tissues_n23 --wo_output 2_RefSeq_Ensembl_Search
conda deactivate
```

</dd></dl>
</p>
</details>
</dd></dl>

<!-- 4. Gene prediction -->
### 4. Gene prediction ###
<dl><dd>
<!-- a. -->
<details><summary><b>a. Run Blastx</b></summary>
<p>
<dl><dd>

```
python -u ./src/Run_blastx.py --dataset 2019_tissues_n23 --blastx ./src/ncbi-blast-2.14.0+/bin/blastx
```

</dd></dl>
</p>
</details>
<!-- b. -->
<details><summary><b>b. Collect hints and run Augustus</b></summary>
<p>
<dl><dd>

```
python -u ./src/Run_Augustus.py --dataset 2019_tissues_n23 --augustus ./src/augustus-3.3.3/bin/augustus
```

</dd></dl>
</p>
</details>
</dd></dl>

<!-- 5. Parse Augustus output -->
### 5. Parse Augustus output ###
<dl><dd>
<!-- a. -->
<details><summary><b>a. Parse Augustus output and run BlastP</b></summary>
<p>
<dl><dd>

```
python -u ./src/Parse_Augustus_and_run_blastp.py --dataset 2019_tissues_n23 --blastp ./src/ncbi-blast-2.14.0+/bin/blastp
```

</dd></dl>
</p>
</details>

<!-- b. -->
<details><summary><b>b. Parse BlastP output</b></summary>
<p>
<dl><dd>

```
python -u ./src/Parse_blastp_XML_output.py --dataset 2019_tissues_n23 --email email@email.com
```

</dd></dl>
</p>
</details> 

<!-- c. Write evidence tables -->
<details><summary><b>c. Write evidence tables</b></summary>
<p>
<dl><dd>

```
python -u ./src/Write_evidence_tables.py --dataset 2019_tissues_n23
```

</dd></dl>
</p>
</details> 
</dd></dl>

<!-- 6. Run MSGF+ against RefSeq+Ensembl+Augustus database --> 
### 6. Run MSGF+ against RefSeq+Ensembl+Augustus database ###
<dl><dd>

>**Note** 
>Requires manual inspection of Augustus predictions and construction of RefSeq+Ensembl+Augustus database.

<!-- a. -->
<details><summary><b>a. Salmon Proteogenomics Project | MSGFPlus_PeptideMatchCMD <a href="https://proteomics2.ucsd.edu/ProteoSAFe/index.jsp">workflow</a></b></summary>
<p>
<dl><dd>

![image](https://drive.google.com/uc?export=view&id=10Cck2wGT6SSD24GkeDlX5yraHKWboKsb) 

</dd></dl>
</p>
</details>

<!-- b. -->
<details><summary><b>b. Download output (<a href="https://proteomics2.ucsd.edu/ProteoSAFe/status.jsp?task=155866fa4e824b7d8dade8221ff4a68e">Example: 2019_tissues_n23 dataset</a>)</b></summary>
<p>
<dl><dd>

![image](https://drive.google.com/uc?export=view&id=1w15YHhLRxq7SlvB9yKWs-jK64vic5cSw)

</dd></dl>
</p>
</details>

<!-- c. -->
<details><summary><b>c. Placement of "ExactMatch_output" & "MSGF_combinedtsv" directories</b></summary>
<p>
<dl><dd>

>**Output:**<br />
>./data/**[dataset]**/5_RefSeq_Ensembl_Augustus_Search/workflow_output/**ExactMatch_output**<br />
>./data/**[dataset]**/5_RefSeq_Ensembl_Augustus_Search/workflow_output/**MSGF_combinedtsv**<br />

* **[dataset]** = 2019_tissues_n23

</dd></dl>
</p>
</details> 

<!-- d. -->
<details><summary><b>d. Compute FDR</b></summary>
<p>
<dl><dd>

```
conda activate py2
python -u ./src/Run_compute_FDR_customDB.py --dataset 2019_tissues_n23 --wo_output 5_RefSeq_Ensembl_Augustus_Search
conda deactivate
```

</dd></dl>
</p>
</details>

<!-- e. -->
<details><summary><b>e. Create normalized PSK expression matrix</b></summary>
<p>
<dl><dd>

```
python -u ./src/Create_PSK_matrix.py --dataset 2019_tissues_n23
```

</dd></dl>
</p>
</details>
</dd></dl>
