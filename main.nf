#!/usr/bin/env nextflow
params.droso="/app/genome_drosophila/genome"
params.myco="/app/genome_mycoplasma/myco_hynorhinis"
params.index="/app/genome_human/human_genome38" //Aqui se indica el path y ep prefijo de los files del genome build.$baseDir es el directorio desde donde se llama a nextflow
 //Aqui especificamos los archivos fastq del directorio
params.outdir="/home/aperonalope/circulares/nextf/results/"  //Aqui especificamos el output directory donde guardaremos algunos files
params.fastqs = "$baseDir/*.fastq.gz"  //Aqui especificamos los archivos fastq del directorio
   
process alineamiento_humano_end_to_end{
	label "HighMemory"
	executor "slurm"
	
	tag "Alineamiento de ${fastq} con el genoma humano"
	input:
	path fastq
	val genoma_humano

	output:
	path "*.sam" 
	
	script:
	"""
	bowtie2 -p 2 -k 1 --no-unal -x ${genoma_humano} -U $fastq -S \$(basename "$fastq" .fastq.gz).sam --score-min C,-16 --mp 2,2 --rdg 0,2 --rfg 0,2
	"""
}

process alineamiento_droso_end_to_end{
	label "HighMemory"
	tag "Alineamiento de ${fastq} con el genoma de drosophila"
	input:
	path fastq
	val genoma_droso

	output:
	path "*.sam" 
	
	script:
	"""
	bowtie2 --sensitive -p 2 -k 1 --no-unal -x ${genoma_droso} -U $fastq -S ${fastq.baseName}.sam 
	"""
}
process alineamiento_myco_end_to_end{
	label "HighMemory"
	tag "Alineamiento de ${fastq} con el genoma de mycoplsama hyorhinis"
	input:
	path fastq
	val genoma_myco

	output:
	path "*.sam" 	
	
	script:
	"""
	bowtie2 --sensitive -p 2 -k 1 --no-unal -x ${genoma_myco} -U $fastq -S ${fastq.baseName}.sam 
	"""
}


process alineamiento_humano_soft{
	label "HighMemory"
	
	tag "Alineamiento de ${fastq} con el genoma humano"
	
	
	input:
	path fastq
	val genoma_humano

	output:
	path "*.csv" 
	
	script:
	"""
	bowtie2 --sensitive -p 2 --score-min C,0 --no-unal -x ${genoma_humano} -U $fastq -S ${fastq.baseName}.sam
	awk 'BEGIN{FS="\\t"} !/^@/{print \$1 " " \$4 " " \$6" "\$3" "\$10}' ${fastq.baseName}.sam |tr "_" " "|sort -k2,2 -k 4,4n > ${fastq.baseName}_limpio.csv
	
	"""
}


process alineamiento_soft_final{
	label "VeryHighMemory"
	publishDir params.outdir , mode: 'copy'
	input:
	path fastq
	val genoma_humano
	
	output:
	path "*.csv"
	
	script:
	"""
	bowtie2 --sensitive-local -p 2 -k 200 --score-min G,2,8 --no-unal -x ${genoma_humano} -U $fastq -S ${fastq.baseName}.sam 
	awk 'BEGIN{FS="\\t"} !/^@/{print \$1 " " \$4 " " \$6" "\$3" "\$10}' ${fastq.baseName}.sam |tr "_" " "|sort -k2,2 -k 4,4n > ${fastq.baseName}_limpio.csv
	
	"""


}

process sam_divide{
	label "LowMemory"
	
	tag "División de $sa"
	input:
	file sa
	
	output:
	file "*.sam" 
	
	script:
	"""
	python3 /app/divide_SAM.py -i $sa

	"""
}


process remove_sequences {
	label "LowMemory"
	tag "Eliminacion de secuencias en $fastq que esten en $sam"

	input:
	tuple val(base), path(fastq), path(sam)
		
	output:
	path "*fastq.gz"
	
	script:
	"""
	python3 /app/remove_duplicated_seq.py -f $fastq -i $sam 

	"""

}

process remove_sequences_myco {
	label "LowMemory"
	tag "Eliminacion de secuencias en $fastq que esten en $sam"

	input:
	tuple val(base), path(fastq), path(sam)
		
	output:
	path "*fastq.gz"
	
	script:
	"""
	python3 /app/remove_duplicated_seq.py -f $fastq -i $sam

	"""

}


process remove_sequences_droso {
	label "LowMemory"
	tag "Eliminacion de secuencias en $fastq que esten en $sam"

	input:
	tuple val(base), path(fastq), path(sam)
		
	output:
	path "*fastq.gz"
	
	script:
	"""
	python3 /app/remove_duplicated_seq.py -f $fastq -i $sam

	"""

}


process repetitive_units{
	label "LowMemory"
	tag "Aislamiento de secuencias repetitivas en ${fastq_gz}"
	input:
	path fastq_gz

	output:
	file "*.fastq.gz"

	script:
	"""
	python3 /app/cirseq_adapt.py -i $fastq_gz -s "r"

	"""
}

process unique_rnas{
	label "LowMemory"
	tag "Aislamiento de especies moleculares únicas en ${fastq_gz}"
	input:
	file fastq_gz

	output:
	path "*"

	script:
	"""
	python3 /app/filter_duplicated.py -i $fastq_gz
	"""
}


process eliminate_rubish{
	label "LowMemory"

	input:
	tuple val(base), path(fastq), path(sam)
	
	output:
	path "*fastq.gz"
	
	script:
	"""
	python3 /app/remove_duplicated_seq_not.py -f $fastq -i $sam
	"""


}
process shifting{
	label "LowMemory"
	
	tag "Obtencion de todas las permutaciones en ${fastq_gz}"
	input:
	file fastq_gz

	output:
	path "*fastq.gz"

	script:
	"""
	python3 /app/shifting.py -i $fastq_gz -s "s"

	"""


}

process shifting2{
	label "LowMemory"
	
	tag "Obtencion de todas las permutaciones en ${fastq_gz}"
	input:
	file fastq_gz

	output:
	path "*fastq.gz"

	script:
	"""
	python3 /app/shifting2.py -i $fastq_gz -s "s"

	"""


}


process deteccion_de_circulares{
	label "LowMemory"
	publishDir params.outdir , mode: 'copy'
	tag "Deteccion de circulares en ${csv}"
	input:
	file csv

	output:
	path "*csv"

	script:
	"""
	python3 /app/circulares.py -i $csv -a /app/anotaciones_homo_sapiens.csv

	"""
}

process limpieza{
	label "LowMemory"
	publishDir params.outdir , mode: 'copy'
	tag "limpieza de ${csv}"
	input:
	file csv
	
	output:
	path "*csv"
	
	script:
	"""
	python3 /app/Cleaning.py -i $csv
	"""


}

workflow {

	Channel.fromPath(params.fastqs).set{gzchannel}
	alineamiento_humano_end_to_end(gzchannel,params.index).set{sam_discard}
	
	gzchannel.map { file -> tuple(file.baseName.replaceAll(/\.fastq$/, ''), file) }.set{gzchannel}
	sam_discard.map{ file -> tuple(file.baseName, file) }.set{sam_discard}
	sam_fastq_join= gzchannel.join(sam_discard)

	

	
	remove_sequences(sam_fastq_join).set{fastqfiltered}
	
	alineamiento_myco_end_to_end(fastqfiltered,params.myco).set{sam_myco}	
	
	fastqfiltered.map{ file -> tuple(file.baseName.replaceAll(/_rm_dupl_filtered\.fastq$/, ''), file)}.set{fastqfiltered2}
	sam_myco.map{ file -> tuple(file.baseName.replaceAll(/_rm_dupl_filtered\.fastq$/, ''),file)}set{sam_discard_myco}

	sam_fastq_join= fastqfiltered2.join(sam_discard_myco)
	remove_sequences_myco(sam_fastq_join).set{fastqfiltered_myco}


	
	
	alineamiento_droso_end_to_end(fastqfiltered_myco,params.droso).set{sam_droso}

	

	fastqfiltered_myco.map{ file -> tuple(file.baseName.replaceAll(/_rm_dupl_filtered_rm_dupl_filtered\.fastq$/, ''), file)}.set{fastqfiltered3}
	sam_droso.map{ file -> tuple(file.baseName.replaceAll(/_rm_dupl_filtered_rm_dupl_filtered\.fastq$/, ''),file)}set{sam_discard_droso}

	sam_fastq_join_final=fastqfiltered3.join(sam_discard_droso)

	remove_sequences_droso(sam_fastq_join_final).set{fastqfiltered_final}

	repetitive_units(fastqfiltered_final).set{consensus}

	consensus
		.flatMap { files ->
			files.findAll { file -> !file.name.endsWith('rotated.fastq.gz') }
		}
		.set { consensus_out }

	unique_rnas(consensus_out).set{consensus_unique}
	consensus_unique
		.flatMap { files ->
			files.findAll { file -> file.name.endsWith('final.fastq.gz') }
		}
		.set {consensus_unique_filtered }

	shifting(consensus_unique_filtered).set{shifted}
	alineamiento_humano_soft(shifted,params.index).set{sam_soft}
	
	deteccion_de_circulares(sam_soft).set{circulares}
	
	limpieza(circulares)

	

	//eyJ0aWQiOiA5MTA3fS4wYzc5MTI0MjNhNmNhMjczZjVjODcyMjlhYzE5OTg1M2NhODM4NTA5
	// for future grafics:https://nf-co.re/docs/contributing/design_guidelines
}
