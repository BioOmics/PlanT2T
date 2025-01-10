#!/usr/bin/bash

spe=$1

if [[ ! -f ${spe}_gene.csv || ! -f ${spe}_transbase.csv || ! -f ${spe}_transfunc.csv ]]; then
    echo "Lack required files...Pleas check"
    exit 255
fi

echo "use plant2t_user;" > ${spe}.sql

## gene --> transbase --> transfunc
echo "
create table \`${spe}_gene\` (
  gene_ID char(100) not null,
  source_ID char(100) not null,
  location char(50) not null,
  start int not null,
  end int not null,
  strand char(20) not null,
  tf_type char(10) null,
  tf_family char(20) null,
  function text null,
  go text null,
  kegg text null,
  primary key (gene_ID)
) ENGINE=MyISAM;
load data local infile '${spe}_gene.csv' into table \`${spe}_gene\` fields terminated by ',' enclosed by '\"' ignore 1 lines;
create table \`${spe}_transfunc\` (
  record_ID int not null auto_increment,
  transcript_ID char(100) not null,
  gene_ID char(100) not null,
  source char(50) not null,
  start int not null,
  end int not null,
  name char(50) not null,
  signature_desc text null,
  Dbxref char(50) null,
  tfend int not null,
  primary key (record_ID)
) ENGINE=MyISAM;
load data local infile '${spe}_transfunc.csv' into table \`${spe}_transfunc\` fields terminated by ',' enclosed by '\"' ignore 1 lines;
create table \`${spe}_transbase\` (
  transcript_ID char(100) not null,
  gene_ID char(100) not null,
  genomic_len int not null,
  transcript_len int not null,
  cds_len int not null,
  pep_len int not null,
  transcript_location char(50) not null,
  gene_strand char(20) not null,
  gene_sequence text not null,
  pep_sequence text not null,
  primary key (transcript_ID)
) ENGINE=MyISAM;
load data local infile '${spe}_transbase.csv' into table \`${spe}_transbase\` fields terminated by ',' enclosed by '\"' ignore 1 lines;
create table \`${spe}_protparam\` (
  SequenceName char(100) not null,
  MolecularWeight float not null,
  Aromaticity float not null,
  InstabilityIndex float not null,
  IsoelectricPoint float not null,
  HelixFraction float not null,
  TurnFraction float not null,
  SheetFraction float not null,
  ReducedCysteinesExtinctionCoefficient int not null,
  OxidizedCysteinesExtinctionCoefficient int not null,
  GRAVY float not null,
  AverageFlexibility float not null,
  ChargeAtpH7_0 float not null,
  primary key (SequenceName)
) ENGINE=MyISAM;
load data local infile '${spe}.protparam.txt' 
into table \`${spe}_protparam\` 
fields terminated by '\\t' 
lines terminated by '\\n' 
ignore 1 lines 
(SequenceName, MolecularWeight, Aromaticity, InstabilityIndex, IsoelectricPoint, HelixFraction, TurnFraction, SheetFraction, ReducedCysteinesExtinctionCoefficient, OxidizedCysteinesExtinctionCoefficient, GRAVY, AverageFlexibility, ChargeAtpH7_0);
" >> ${spe}.sql

while IFS=$'\t' read -r db_ID db_version id name common_name full_name img_src collected_date Ploidy total_length fragments_num N50 N_count Gaps T2TCount T2NCount N2NCount gene_num transcript_num busco TR TF GC EndStop MidStop NoStop ncbi_taxonomy_ID Phylum Class Order Family Genus QV LAI Email Author Organization DOI Message assem_source GenomeSurvey AssemblyByHiFi AssemblyByONT AssemblyChloroplast Polish HiCScaffolding GapFilling QualityAssessment TelomereIdentify CentromereIdentify Subgenome RepeatAnno CodingAnno
do
  echo "INSERT INTO tb_genome (db_ID, db_version, id, name, common_name, full_name, img_src, collected_date, Ploidy, total_length, fragments_num, N50, N_count, Gaps, T2TCount, T2NCount, N2NCount, gene_num, transcript_num, busco, TR, TF, GC, EndStop, MidStop, NoStop, ncbi_taxonomy_ID, Phylum, Class, \`Order\`, Family, Genus, QV, LAI, Email, Author, Organization, DOI, Message, assem_source, GenomeSurvey, AssemblyByHiFi, AssemblyByONT, AssemblyChloroplast, Polish, HiCScaffolding, GapFilling, QualityAssessment, TelomereIdentify, CentromereIdentify, Subgenome, RepeatAnno, CodingAnno)
  VALUES ('$db_ID', '$db_version', '$id', '$name', '$common_name', '$full_name', '$img_src', '$collected_date', '$Ploidy','$total_length', $fragments_num, $N50, $N_count, $Gaps, $T2TCount, $T2NCount, $N2NCount, $gene_num, $transcript_num, '$busco', $TR, $TF, '$GC', $EndStop, $MidStop, $NoStop, '$ncbi_taxonomy_ID', '$Phylum', '$Class', '$Order', '$Family', '$Genus', '$QV', '$LAI', '$Email', '$Author', '$Organization', '$DOI', '$Message', '$assem_source', '$GenomeSurvey', '$AssemblyByHiFi', '$AssemblyByONT', '$AssemblyChloroplast', '$Polish', '$HiCScaffolding', '$GapFilling', '$QualityAssessment', '$TelomereIdentify', '$CentromereIdentify', '$Subgenome', '$RepeatAnno', '$CodingAnno');" >> ${spe}.sql
done < genome.summary.txt

cat ${spe}.sql | mysql -u chendijun -pchen2019@525!
