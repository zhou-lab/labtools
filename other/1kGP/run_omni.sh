export SV_DIR=/scratch/bcb/wzhou1/tools/genomestrip/svtoolkit

# Omni chip

java -Xmx4g -cp /scratch/bcb/wzhou1/tools/genomestrip/svtoolkit/lib/SVToolkit.jar:/scratch/bcb/wzhou1/1kGP/2013_11_IRS/gatk/GenomeAnalysisTK-2.7-4-g6f46d11/GenomeAnalysisTK.jar:/scratch/bcb/wzhou1/tools/genomestrip/svtoolkit/lib/gatk/Queue.jar org.broadinstitute.sv.main.SVAnnotator -A IntensityRankSum -R /scratch/bcb/wzhou1/reference/human_g1k_v37.fasta -vcf $1 -arrayIntensityFile /scratch/bcb/wzhou1/1kGP/2013_11_IRS/ALL.genome.Omni25_probe_intensity_matrix_2141samples.20120131.dat -irsUseGenotypes true -writeReport TRUE -reportFile $1.omni.dat
