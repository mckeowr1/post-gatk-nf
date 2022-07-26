/*
==================================
~ > *                        * < ~
~ ~ > *                    * < ~ ~
~ ~ ~ > *  ANNOTATE VCF  * < ~ ~ ~
~ ~ > *                    * < ~ ~
~ > *                        * < ~
==================================
*/


/*
------------ Extract ancestor strain from the VCF and make bed file for annotations 
*/

process extract_ancestor_bed {

    label 'postgatk'

    publishDir "${params.output}/ANNOTATE_VCF", mode: 'copy'

    cpus 1

    input:
      tuple file(vcf), file(vcfindex)

    output:
      tuple file("ANC.bed.gz"), file("ANC.bed.gz.tbi")

      """
        bcftools query --samples ${params.anc} -f '%CHROM\\t%POS\\t%END\\t[%TGT]\\n' ${vcf} |\\
        awk -F"/" '\$1=\$1' OFS="\\t" |\\
        awk '{print \$1, \$2 = \$2 - 1, \$3, \$4}' OFS="\\t" |\\
        bgzip > ANC.bed.gz

        tabix ANC.bed.gz
        echo "ANCESTOR DONE"
      """
}

/*
------------ Annotate small variant VCF  
*/

process annotate_small_vcf {

    publishDir "${params.output}/ANNOTATE_VCF", mode: 'copy'

    label 'postgatk'

    // conda '/projects/b1059/software/conda_envs/vcffixup'

    cpus 1

    input:
      tuple file(vcf), file(vcfindex), file("ANC.bed.gz"), file("ANC.bed.gz.tbi"), file(pop) //, val(pop), val(maf), val(sm)

    output:
      tuple file("Ce330_annotated.vcf.gz"), file("Ce330_annotated.vcf.gz.tbi")


      """
        # get vcfanno files
        cp ${workflow.projectDir}/input_files/annotations/${params.species}/* .
        cat ${params.vcfanno_config} | sed 's/species/${params.species}/' > anno_config.toml

        bcftools view -S ${pop} ${vcf} -Oz -o population-filtered.vcf.gz 

        vcfanno anno_config.toml population-filtered.vcf.gz |\\
        awk '\$0 ~ "#" || \$0 !~ "Masked" {print}' |\\
        vcffixup - |\\
        bcftools filter -i N_MISSING=0 -Oz -o Ce330_annotated.vcf.gz

        tabix -p vcf Ce330_annotated.vcf.gz
      """
}

/*
------------ Annotate snv vcf with genetic distance
*/

process annotate_vcf_genetic_distance {

    publishDir "${params.output}/ANNOTATE_VCF", mode: 'copy'

    label 'postgatk'
    memory 20.GB

    // conda '/projects/b1059/software/conda_envs/vcffixup'

    cpus 1

    input:
      tuple file(vcf), file(vcfindex)//, val(pop), val(maf), val(sm)

    output:
      tuple file("Ce330_annotated.vcf.gz"), file("Ce330_annotated.vcf.gz.tbi")


      """
        # get vcfanno files
        cp ${workflow.projectDir}/input_files/annotations/${params.species}/* .
        cat ${params.vcfanno_config} | sed 's/species/${params.species}/' > anno_config.toml


        vcfanno anno_config.toml ${vcf} |\\
        awk '\$0 ~ "#" || \$0 !~ "Masked" {print}' |\\
        vcffixup - |\\
        bcftools filter -i N_MISSING=0 -Oz -o Ce330_annotated.vcf.gz

        tabix -p vcf Ce330_annotated.vcf.gz
      """
}

/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  Run PCA and DAPC  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/

/*
------------ Mask Divergent Positions
*/

process mask_divergent{
  tag{"MASK DIVERGENT SITE"}
  
  label 'postgatk'
  
  memory 20.GB
  cpus 4
  
  publishDir "${params.output}/DIVMASK", mode: 'copy'
  
  input:
    tuple file(vcf), file(vcfindex), file(regions)
  
  
  output:
    tuple file("masked.vcf.gz"), file("masked.vcf.gz.tbi")
  
  """
  bcftools view --threads 3 -O z -T ${regions} ${vcf} > masked.vcf.gz  
  
  bcftools index --threads 3 --tbi masked.vcf.gz 
  
  """
  
}

process filter_markers {

  tag {"filter_markers"}

  label 'postgatk'
  
  memory 30.GB
  cpus 6
  time '3h'

  // conda '/projects/b1059/software/conda_envs/vcffixup'

  // publishDir "${params.output}/EIGESTRAT/LD_${test_ld}/INPUTFILES", mode: 'copy'

  input:
    tuple file(vcf), file(vcfindex), val("test_ld")

  output:
    tuple val(chrom), val(test_ld), file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind"), file("plink.prune.in"), \
    file ("markers.txt"), file ("sorted_samples.txt"), file ("PCA.vcf.gz"), file ("PCA.vcf.gz.tbi")


    """

    bcftools view -Ou --regions ${chrom} ${vcf} |\\
    bcftools norm --threads 5 -m + -Oz -o ce_norm.vcf.gz

    tabix -p vcf ce_norm.vcf.gz

    plink --vcf ce_norm.vcf.gz --biallelic-only --set-missing-var-ids @:# --indep-pairwise 50 10 ${test_ld} --allow-extra-chr 

    plink --vcf ce_norm.vcf.gz --biallelic-only --set-missing-var-ids @:# --extract plink.prune.in --geno --recode12 --out eigenstrat_input --allow-extra-chr

    awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
    sort -k1,1d -k2,2n > markers.txt

    bcftools query -l ce_norm.vcf.gz |\\
    sort > sorted_samples.txt 

    bcftools view -S sorted_samples.txt -R markers.txt ce_norm.vcf.gz -Oz -o PCA.vcf.gz
    
    tabix -p vcf PCA.vcf.gz

    bcftools view -S sorted_samples.txt -R markers.txt ce_norm.vcf.gz |\\
    bcftools query -f '%CHROM\\t%CHROM:%POS\\t%POS\\t%REF\\t%ALT\\n' |\\
    sed 's/^III/3/g' |\\
    sed 's/^II/2/g' |\\
    sed 's/^IV/4/g' |\\
    sed 's/^I/1/g' |\\
    sed 's/^V/5/g' > eigenstrat_input.pedsnp      

    cut -f-6 -d' ' eigenstrat_input.ped |\\
    awk '{print 1, \$2, \$3, \$3, \$5, 1}'  > eigenstrat_input.pedind

    echo "rerun"
    """

}

/*
------------ Prepare files for EIGENSTRAT
*/

process vcf_to_eigstrat_files {

  tag {"PREPARE EIGENSTRAT FILES"}

  label 'postgatk'
  
  memory 30.GB
  cpus 6
  time '3h'

  // conda '/projects/b1059/software/conda_envs/vcffixup'

  publishDir "${params.output}/EIGESTRAT/LD_${test_ld}/INPUTFILES", mode: 'copy'

  input:
    tuple file(vcf), file(vcfindex), val(test_ld)

  output:
    tuple file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind"), file("plink.prune.in"), \
    file ("markers.txt"), file ("sorted_samples.txt"), file ("PCA.vcf.gz"), file ("PCA.vcf.gz.tbi")


    """

    bcftools view -Ou --regions I,II,III,IV,V,X ${vcf} |\\
    bcftools norm --threads 5 -m + -Oz -o ce_norm.vcf.gz

    tabix -p vcf ce_norm.vcf.gz

    plink --vcf ce_norm.vcf.gz --biallelic-only --set-missing-var-ids @:# --indep-pairwise 50 10 ${test_ld} --allow-extra-chr 

    plink --vcf ce_norm.vcf.gz --biallelic-only --set-missing-var-ids @:# --extract plink.prune.in --geno --recode12 --out eigenstrat_input --allow-extra-chr

    awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
    sort -k1,1d -k2,2n > markers.txt

    bcftools query -l ce_norm.vcf.gz |\\
    sort > sorted_samples.txt 

    bcftools view -S sorted_samples.txt -R markers.txt ce_norm.vcf.gz -Oz -o PCA.vcf.gz
    
    tabix -p vcf PCA.vcf.gz

    bcftools view -S sorted_samples.txt -R markers.txt ce_norm.vcf.gz |\\
    bcftools query -f '%CHROM\\t%CHROM:%POS\\t%POS\\t%REF\\t%ALT\\n' |\\
    sed 's/^III/3/g' |\\
    sed 's/^II/2/g' |\\
    sed 's/^IV/4/g' |\\
    sed 's/^I/1/g' |\\
    sed 's/^V/5/g' > eigenstrat_input.pedsnp      

    cut -f-6 -d' ' eigenstrat_input.ped |\\
    awk '{print 1, \$2, \$3, \$3, \$5, 1}'  > eigenstrat_input.pedind

    echo "rerun"
    """

}


/*
------------ Get a list of Markers that pass filtering thresholds 
*/

process get_passing_variants{
  tag {"PREPARE EIGENSTRAT FILES"}

  label 'postgatk'
  publishDir "${params.output}/EIGESTRAT/LD_${test_ld}/PASSING/", mode: 'copy'
  
  input:
    tuple file(vcf), file(vcfindex), val("test_ld")
  
  output:
    tuple val(test_ld), file("ce_norm.vcf.gz"), file ("ce_norm.vcf.gz.tbi"), file ("markers.txt"), file ("sorted_samples.txt")

  """
    bcftools view -Ou --regions I,II,III,IV,V,X ${vcf} |\\
    bcftools norm --threads 5 -m + -Oz -o ce_norm.vcf.gz

    tabix -p vcf ce_norm.vcf.gz

    plink --vcf ce_norm.vcf.gz --biallelic-only --set-missing-var-ids @:# --indep-pairwise 50 10 ${test_ld} --allow-extra-chr 

    plink --vcf ce_norm.vcf.gz --biallelic-only --set-missing-var-ids @:# --extract plink.prune.in --geno --recode12 --out eigenstrat_input --allow-extra-chr

    awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
    sort -k1,1d -k2,2n > markers.txt

    bcftools query -l ce_norm.vcf.gz |\\
    sort > sorted_samples.txt 

  """
}


/*
------------ Filter VCF to markers 
*/


process filter_vcf{
  publishDir "${params.output}/EIGESTRAT/LD_${test_ld}/VCFs/", mode: 'copy'

  tag {"PREPARE EIGENSTRAT FILES"}

  label 'postgatk'
  input: 
    tuple val("chrom"), val("test_ld"), file(marker_vcf), file(index), file(markers_list), file(strains_list)
  
  output:
    tuple val(test_ld), path("${chrom}_${test_ld}_filtered_vcf.gz"), file("${chrom}_${test_ld}_filtered_vcf.gz")
  
  """
  bcftools view --regions ${chrom} -Oz -o ${chrom}_vcf.gz ${marker_vcf} 
  tabix -p vcf ${chrom}_vcf.gz
  bcftools view -S ${strains_list} -R ${markers_list} -Oz -o ${chrom}_${test_ld}_filtered_vcf.gz ${chrom}_vcf.gz
  tabix -p vcf ${chrom}_${test_ld}_filtered_vcf.gz
  """
}

/*
------------ Concat Chrom VCFS 
*/


process concat_vcf{
  input: 
    tuple val(test_ld), file(marker_vcf), file(index), file(markers_list), file(strains_list)
  
  output:
    tuple val(test_ld), file(pca_chrom_vcf), file(pca_chrom_index)
  
  """
  bcftools view -Ou --regions ${chrom} ${marker_vcf} |\\
  bcftools view -S ${strains_list} -R ${markers_list} -Oz -o ${chrom}_${test_ld}_filtered_vcf.gz 
  tabix -p vcf ${chrom}_${test_ld}_filtered_vcf.gz
  """
}


/*
------------ Run EIGENSTRAT without removing outlier strains
*/

process run_eigenstrat_no_outlier_removal {

  publishDir "${params.output}/EIGESTRAT/LD_${test_ld}/NO_REMOVAL/", mode: 'copy'

  label 'pca'
  
  memory 50.GB

  // conda '/projects/b1059/software/conda_envs/vcffixup'

  input:
    tuple val("ld"), file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind"), file("plink.prune.in"), \
    file ("markers.txt"), file ("sorted_samples.txt"), file ("PCA.vcf.gz"), file ("PCA.vcf.gz.tbi"), file(eigenparameters)

  output:
    tuple val(ld), file("eigenstrat_no_removal.evac"), file("eigenstrat_no_removal.eval"), file("logfile_no_removal.txt"), \
    file("eigenstrat_no_removal_relatedness"), file("eigenstrat_no_removal_relatedness.id"), file("TracyWidom_statistics_no_removal.tsv")


    """

    smartpca -p ${eigenparameters} > logfile_no_removal.txt

    sed -n -e '/Tracy/,\$p' logfile_no_removal.txt |\
    sed -e '/kurt/,\$d' |\
    awk '\$0 !~ "##" && \$0 !~ "#" {print}' |\
    sed -e "s/[[:space:]]\\+/ /g" |\
    sed 's/^ //g' |\
    awk 'BEGIN{print "N", "eigenvalue", "difference", "twstat", "p-value", "effect.n"}; {print}' OFS="\\t" |\
    awk -F" " '\$1=\$1' OFS="\\t" > TracyWidom_statistics_no_removal.tsv
    """

}

/*
------------ Run EIGENSTRAT with removing outlier strains
*/

process run_eigenstrat_with_outlier_removal {

  // conda '/projects/b1059/software/conda_envs/vcffixup'

  label 'pca'
  memory 50.GB

  publishDir "${params.output}/EIGESTRAT/LD_${test_ld}/OUTLIER_REMOVAL/", mode: 'copy'

  input:
    tuple val("ld"), file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind"), file("plink.prune.in"), \
    file ("markers.txt"), file ("sorted_samples.txt"), file ("PCA.vcf.gz"), file ("PCA.vcf.gz.tbi"), file(eigenparameters)

  output:
    tuple val(ld), file("eigenstrat_outliers_removed.evac"), file("eigenstrat_outliers_removed.eval"), file("logfile_outlier.txt"), \
    file("eigenstrat_outliers_removed_relatedness"), file("eigenstrat_outliers_removed_relatedness.id"), file("TracyWidom_statistics_outlier_removal.tsv")

   
    """
    smartpca -p ${eigenparameters} > logfile_outlier.txt

    sed -n -e '/Tracy/,\$p' logfile_outlier.txt |\
    sed -e '/kurt/,\$d' |\
    awk '\$0 !~ "##" && \$0 !~ "#" {print}' |\
    sed -e "s/[[:space:]]\\+/ /g" |\
    sed 's/^ //g' |\
    awk 'BEGIN{print "N", "eigenvalue", "difference", "twstat", "p-value", "effect.n"}; {print}' OFS="\\t" |\
    awk -F" " '\$1=\$1' OFS="\\t" > TracyWidom_statistics_outlier_removal.tsv
    """

}


/*
------------ Run HTML report for PCA analysis
*/

process HTML_report_PCA {

  label 'R'

  // conda '/projects/b1059/software/conda_envs/cegwas2-nf_env'

  publishDir "${params.output}/", mode: 'copy'


  input:
  tuple val("ld"), file("eigenstrat_no_removal.evac"), file("eigenstrat_no_removal.eval"), file("logfile_no_removal.txt"), \
    file("eigenstrat_no_removal_relatedness"), file("eigenstrat_no_removal_relatedness.id"), file("TracyWidom_statistics_no_removal.tsv"), \
    file("eigenstrat_outliers_removed.evac"), file("eigenstrat_outliers_removed.eval"), file("logfile_outlier.txt"), \
    file("eigenstrat_outliers_removed_relatedness"), file("eigenstrat_outliers_removed_relatedness.id"), \
    file("TracyWidom_statistics_outlier_removal.tsv"), file(pca_report), file(pca_template)


  output:
   tuple file("pca_report.Rmd"), file('pca_template.Rmd'), file("*.html")


  """
  # prepare for report
  cat ${pca_report} | \\
  sed "s+LD_VALUE+${ld}+" | \\
  sed "s+EIGESTRAT/{ld}/NO_REMOVAL/++g" | \\
  sed "s+EIGESTRAT/{ld}/OUTLIER_REMOVAL/++g" > pca_report_LD_${ld}.Rmd

  Rscript -e "rmarkdown::render('pca_report_LD_${ld}.Rmd')"
  """


}

//  cat pca_report.Rmd | sed 's+pca_template.Rmd+${workflow.projectDir}/bin/pca_template.Rmd+' > new_pca_report.Rmd



