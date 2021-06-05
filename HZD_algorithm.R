#VarGenius-HZD searches rare homozygous and hemizygous deletions in targeted sequencing
#Copyright (C) 2021 Francesco Musacchia

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


#Description: This R script performs the detection of relevant Homozygous/Hemizygous deletions from 
#results of coverage analysis

#Input:
# 1. a table with: samplename\tinput_folder\n must be inputed that is needed to know where are the results of coverage for each sample
# 2. the enrichment kit name (full name without .bed)
# 3. the suffix of the exonic intervals used (UCSC_genes)
# 4. the suffix to use for the non covered regions (noncov)
# 5. the name to use for the final output


#ADDS ON:
# 1.
# After the creation of the suspect table, outliers should be removed by simply deciding a threshold.
# For example all those that have the double of CNVs compared with the next one in the order per CNVs number

# 2.
# Removal of those deleted exons that actually are at the end of the splicing and could be removed for some reason


#Gets the parameters in input
args <- commandArgs(trailingOnly = TRUE)

length(args)
if (length(args) < 1){
	stop("USAGE: R CMD BATCH --no-save --no-restore '--args <samples_tab_f> <bamfile_suff> <enrichmentkit_name> <exons_suff> <noncov_suff> <final_out_ext>' CNV_analysis.R")
}

#This table should be: 
#samplename\tpath
samples_tab_f = args[1]
samples_tab_f

#target limited to exons name without extension
targetexons_name = args[2]
targetexons_name

#The suffix for the genomecov output tables
exons_suff = args[3]
exons_suff

#The suffix for the non coovered exons tables
noncov_suff = args[4]
noncov_suff

#output folder
out_f = args[5]
out_f

#output suffix
out_suff = args[6]
out_suff


nc_suffix = paste(targetexons_name,noncov_suff,sep="_")
#gen_cov_suffix =  paste(paste(enrichmentkit_name,exons_suff,sep="_"),"bed",sep=".")
gen_cov_suffix =  paste(targetexons_name,"bed",sep=".")

#General parameters
#Type of analysis and CNV type to detect. As far as 2020 this tool can detect only DELetions
analtype = "single"#"family"
cnv_type_suff = "DEL" #These are deletions
#I use here 2 because I consider that a maximum of 2 samples have the same non-covered region
max_samples_with_lowcov = 2


#Open the samples table with included paths for coverage analysis output
samples_tab <- read.delim(samples_tab_f,stringsAsFactors=F,head=T,quote="",sep="\t")



#Splits the string using the underscore (_) and removes the last element.
#Then concatenates again the parts. That is the familyname!
get_family_name = function (samplename) {
	#Check if the name contains _P meaning that this is a proband
	if(grepl('_P',samplename)){
		list = unlist(strsplit(samplename,"_"))
		a<-list[1:length(list)-1]
		familyname <- paste(a,collapse='_')
	}else{
		familyname <- samplename
	}
	familyname
}


#####################
####STEP 1 #### COLLECTION OF ALL NON-COVERED REGIONS AND CORRESPONDING GENE ANNOTATION 
#####################



#This first part of the script has the goal to read all the results of the 
#non covered regions extraction operated using bedtools
#Family analysis: the analysis can be conducted using the data of the parents, if available
#These can be used for the comparison of the coverages per-exon
interval_2_gene <- data.frame(compid=character(),
                 Gene=character(),
                 stringsAsFactors=FALSE) 
                 
nc_genes_all <- data.frame(Date=as.Date(character()),
                 File=character(), 
                 User=character(), 
                 stringsAsFactors=FALSE) 

samplenum = 1
for ( samples_count in 1:nrow(samples_tab) ){
	
	samplename <- samples_tab[samples_count,c("samplename")]
	samplename
	bampath <- samples_tab[samples_count,c("path")]
	input_fold<-dirname(bampath)	
	
	nc_genes_t <- paste(input_fold,paste(samplename,nc_suffix,sep="_"),sep="/")
	print (nc_genes_t)
	
	#Check if file exists and is not empty
	if (file.info(nc_genes_t)$size > 0  & file.exists(nc_genes_t) ){
		nc_genes <- read.delim(nc_genes_t,stringsAsFactors=F,head=F,quote="",sep="\t")
		
		#When the decimal is separated with comma, change it
		nc_genes$V7<-sub(",",".",nc_genes$V7)
		
		#Get a composite id with chr_start_end of the exon
		nc_genes$compid <-  paste(nc_genes$V1,paste(nc_genes$V2,nc_genes$V3,sep="-"),sep=":")
		nc_genes_filt <- nc_genes[,c("compid","V7")]
		
		#Set column names
		colnames(nc_genes_filt) <- c("exon",samplename)
		
		if (nrow(nc_genes_filt) > 0){

			#Put it as numeric
			nc_genes_filt[[samplename]] <- as.numeric(nc_genes_filt[[samplename]])
			#Define the complete file if it is empty
			if (samplenum == 1){
				nc_genes_all <- nc_genes_filt
				interval_2_gene <- nc_genes[,c("V11","compid")]
				#colnames(genes_corresp) <- c("Chr","Start","End","Gene")
			}else{
				#Add the non covered regions for this sample to the overall table
				nc_genes_all <- merge(nc_genes_all,nc_genes_filt,by="exon",all=T)
				
				#From the file of non covered regions take also the gene annotation correspondingly with 
				#those interval. Since exactly same interval could be present multiple times for different samples
				#then check if the same interval does not exists yet in this set
				#E.g. V11                     compid
				#     EPHB2,MIR4253,MIR4684   chr1:23037475-23037536
				nc_genes_2 <-nc_genes[!(nc_genes$compid %in% interval_2_gene$compid),c("V11","compid")]
				interval_2_gene <- rbind(interval_2_gene,nc_genes_2)
			}
			samplenum = samplenum + 1			
		}		
	}else{
		print(paste(nc_genes_t,"cannot be used",sep =" "))
	}	
}


#####################
####STEP 2 #### GET SUSPECT NON COVERED REGIONS 
#####################


get_suspect_regions = function (suspect_table_f,nc_genes_all_temp) {
	
	
	#This second part evaluates if an interval is suspect 
	#Only if there are few samples (threshold:max_samples_with_lowcov - suggested max 2) 
	#with a certain low coverage, the interval is considered suspect, otherwise is a systematic error or difficulty of sequencing
	#(e.g. high GC content regions)
	suspect_table <- data.frame(Interval=character(),
					 samplename=character(), 
					 Value=character(), 
					 stringsAsFactors=FALSE) 
					 
	#Now count the occurrencies of percentages.  Equal percentages are most probably regions which cannot be
	#seen by the sequencer
	out_name = paste(targetexons_name,"suspect_intervals.txt",sep="_")
	cat(paste("interv_uniq","sam_name_uniq","value",sep="\t"),file=out_name,append = FALSE)
	cat("\n",file=out_name,append = TRUE)


	interv_num = 0
	#Define the complete file if it is empty
	if (nrow(nc_genes_all_temp) > 0){
		#Where it finds NA means that the exon is well covered, hence put 1
		nc_genes_all_temp[is.na(nc_genes_all_temp)] <- 1
		nc_genes_all_t <- t(nc_genes_all_temp)
		#For each interval evaluate the proportion of coverage
		for ( interval in 2:ncol(nc_genes_all_t)){
			#Count the number of time each percentage (means amount of coverage) is present
			freq_table <-as.data.frame(ftable(nc_genes_all_t[,interval]))
			#It counts also the interval name. Hence I remove this count
			freq_table_ok <- freq_table[1:nrow(freq_table)-1,]
			#Leave only frequencies which could be useful
			#This step is really important because it actually filters only those intervals which have
			#a suspect coverage removing the following percentages:
			# - value = 1 : that means the exon is covered for all samples, it has no importance
			# - value = same value more than X times: means that this region is systematically non-covered in several samples
			freq_table_ok_f <- freq_table_ok[as.vector(freq_table_ok$Var1) < 1 & as.vector(freq_table_ok$Freq) < max_samples_with_lowcov,]
			
			#Only if there are very few sample with a low coverage,consider this interval
			#This is an additional control: there may be two different low percentages, but we consider a maximum with the same threshold
			#Hence if max_samples_with_lowcov=2, a maximum of 4 total samples could have low coverage here.
			if (nrow(freq_table_ok_f) <= max_samples_with_lowcov  && nrow(freq_table_ok_f) > 0 ){
				#Get the value
				value <- as.numeric(as.vector(freq_table_ok_f$Var1)[1])
				#Get the name of the samples which has this "unique" low coverage
				sel_samplename <- colnames(nc_genes_all_temp[which(nc_genes_all_temp[interval,]==value)])
				#Get the interval for which this sample has a SV
				interv_uniq <- nc_genes_all_temp$exon[interval]

				#Add interval,samplename,value to a list of suspect SV
				suspect_table[nrow(suspect_table) + 1,] = list(interv_uniq,sel_samplename,value)
			}
		}
	}	
	#Print the suspect table. Useful to check for outliers
	write.table(suspect_table,file=suspect_table_f,sep="\t",row.names=F, quote=F)	
	suspect_table
}


#Here I tried to verify if it was possible to make a loop that until no outlier are detected,
#it continue sto remove them and restart the suspect search. What I observed is that 
#there are different outliers after different rounds of filtering. Finally I decided to
#execute a single and uniq round of filtering
toremove<-list()

nc_genes_all_temp <- nc_genes_all
#Remove outliers and re launch the search of suspect
suspect_table_f <- paste(out_f,paste(targetexons_name,out_suff,"suspect_table.txt",sep="_"),sep="/")
suspect_table <- get_suspect_regions(suspect_table_f,nc_genes_all_temp)


##Open the table and calculate number of CNVs per-sample to filter those samples
##which are more than 3 stdev from the mean
##suspect_table <- read.delim(suspect_table_f,stringsAsFactors=F,head=T,quote="",sep="\t")
#cf <- as.data.frame(ftable(suspect_table[,"samplename"]))
#cf_f <- cf[cf$Freq > (sd(cf$Freq) * 3), ]
#a <- list(as.vector(cf_f$Var1))
#toremove <- append(toremove,a)

##If the number of outliers is greater then one, filter out the outlier and restart the suspect search
#if (nrow(cf_f) > 0){
#	nc_genes_all_f <- nc_genes_all_temp[ , -which(names(nc_genes_all_temp) %in% toremove)]
#	nc_genes_all_temp <- nc_genes_all_f	
#	suspect_table <- get_suspect_regions(suspect_table_f,nc_genes_all_temp)
#}
	
	
#####################
####STEP 3 #### ADD RAW COVERAGE OF OTHER SAMPLES (PREFERABLY PARENTS, too)
#####################

#The suspect non covered exons later should be compared with the other interval:
#those of the parents and also we add the average number of reads present in 
#the other samples, because it can happen that other samples BAM cover the interval but
#with very few reads	#
#if ( nrow(suspect_table) > 0){
	
	#Get the compid (intervals) to be used
	sel_intervals <- as.data.frame(unique(suspect_table$Interval))
	colnames(sel_intervals) <- c("compid")
	geno_cov_all <- data.frame(Date=as.Date(character()),
					 File=character(), 
					 User=character(), 
					 stringsAsFactors=FALSE) 	
							 
	samplenum = 1
	#Get the information for the suspect intervals from the genomeCoverageBed files for both probands and parents
	#If it is a TRIO, the parents will be inserted only once, if >3 then could insert parents data more than once
	#Hence here we need a check: the inserted family will have a vector with already inserted families
	consid_families = c("NONE")
	#for (samplename in samples){
	
	for ( samples_count in 1:nrow(samples_tab) ){
		samplename <- samples_tab[samples_count,c("samplename")]
		bampath <- samples_tab[samples_count,c("path")]
		input_fold<-dirname(bampath)	
				
		print (paste("Getting raw coverage for: ",samplename,sep =" "))
		#Check if the user wants the analysis per-family
		if (analtype == 'family'){				
			family_name = get_family_name(samplename)
			#If the parents have not yet been considered..
			if (!grepl(family_name,consid_families)){
				components = c(samplename,paste(family_name,"M",sep="_"),paste(family_name,"F",sep="_"))
				consid_families[length(consid_families)+1] = family_name					
			}else{components = c(samplename)}
		}else{components = c(samplename)}
		
		#For each component of the family, get the genomeCovBed table 
		for (component in components) {
			geno_cov_f <- paste(input_fold,paste(component,gen_cov_suffix,sep="_"),sep="/")
			#Check that the file exists and is non empty
			if (file.info(geno_cov_f)$size > 0 & file.exists(geno_cov_f)){
				geno_cov <- read.delim(geno_cov_f,stringsAsFactors=F,head=F,quote="",sep="\t")
				#Get a composite id with chr_start_end of the exon
				geno_cov$compid <- paste(geno_cov$V1,paste(geno_cov$V2,geno_cov$V3,sep="-"),sep=":")
				#Filter only the needed intervals
				geno_cov_needed <- merge(geno_cov,sel_intervals,by="compid",all=F)
				#Select columns to investigate
				geno_cov_filt <- geno_cov_needed[,c("compid","V4")]
				colnames(geno_cov_filt) <- c("compid",component)
				#Define the complete file if it is empty
				if (samplenum == 1){
					geno_cov_all <- geno_cov_filt	
				}else{
					#Merge with the all file	
					geno_cov_all <- merge(geno_cov_all,geno_cov_filt,by="compid",all=T)
				}	
				samplenum = samplenum+1						
			}else{
				print(paste(geno_cov_f,"cannot be used",sep =" "))
			}

		}
	}
	
	
	#Convert to numeric
	geno_cov_all_n<-as.matrix(sapply(geno_cov_all[, 2:ncol(geno_cov_all)], as.numeric))
	#re-Assign row names
	rownames(geno_cov_all_n)<-geno_cov_all$compid
	#Get the average for eachrow 
	geno_cov_all_n_st <- cbind(geno_cov_all_n,mean=rowMeans(geno_cov_all_n),compid=geno_cov_all$compid)

	#This last step is needed to print a final table with the suspect
	#inervals and the read coverage from genomecov for the parents if present
	
	#samplename - interval value - prob_reads - fath_reads - moth_reads -avg_reads
	sel_intervals<-as.matrix(sel_intervals)
					
	final_table <- data.frame(
		 samplename=character(), 
		 compid=character(),
		 value=character(),
		  preads=character(),  
		  freads=character(),  
		  mreads=character(),  
		  avgreads=character(),  
		 stringsAsFactors=FALSE)
	 #Go through the intervals to be considered
	 tot_intervals <- length(sel_intervals)
	for (interv_num in 1:tot_intervals ){
		#Get the samplename
		samplenames <- suspect_table[suspect_table$Interval == sel_intervals[interv_num],c("samplename")]
		
		for (samplename in samplenames){
			print(paste(interv_num,"of",tot_intervals,"suspect interval: ",sel_intervals[interv_num]," - Sample:",samplename,sep=" "))
			#Get the family name
			family_name = get_family_name(samplename)
			#Get the value for the interval from the suspect table
			value <- suspect_table[suspect_table$Interval == sel_intervals[interv_num],c("Value")]
			#Get the number of reads for this interval for the proband and the parents
			preads <- geno_cov_all[geno_cov_all$compid == sel_intervals[interv_num],samplename]	
			freads <- geno_cov_all[geno_cov_all$compid == sel_intervals[interv_num],c(paste(family_name,"F",sep="_"))]	
			mreads <- geno_cov_all[geno_cov_all$compid == sel_intervals[interv_num],c(paste(family_name,"M",sep="_"))]	

			avgreads <- geno_cov_all_n_st[as.data.frame(geno_cov_all_n_st)$compid == sel_intervals[interv_num],c("mean")]
			
			#All the information must not null, otherwise I will get an error during writing
			if (length(preads)){
				preads_str <- preads
			}else{
				preads_str <- '-'
			}
			if (length(freads)){
				freads_str <- freads
			}else{
				freads_str <- '-'
			}
			if (length(mreads)){
				mreads_str <- mreads
			}else{
				mreads_str <- '-'
			}			
			if (length(avgreads)){
				avgreads_str <- avgreads
			}else{
				avgreads_str <- '-'
			}	
			#Since this are all deletions I am adding the cnv_type_suff
			outcompid = paste(sel_intervals[interv_num],cnv_type_suff,sep="_")
			final_table[nrow(final_table) + 1,] = list(samplename,outcompid,value,preads_str,mreads_str,freads_str,avgreads_str)				
			
		}

	}		
	
	final_table <- final_table[order(final_table$samplename),]
	
	#Add chr, start, end cnv_type
	final_table$chr <- sapply(strsplit(as.character(final_table$compid),':'), "[", 1)
	temp<-sapply(strsplit(as.character(final_table$compid),':'), "[", 2)
	final_table$start <- as.numeric(sapply(strsplit(as.character(temp),'-'), "[", 1))
	temp2<-sapply(strsplit(as.character(temp),'-'), "[", 2)
	final_table$end <- as.numeric(sapply(strsplit(as.character(temp2),'_'), "[", 1))
	final_table$cnv_type = sapply(strsplit(as.character(temp2),'_'), "[", 2)
	#Print output	
	write.table(final_table,file=paste(out_f,paste(targetexons_name,out_suff,"final_suspect_table.txt",sep="_"),sep="/"),sep="\t",row.names=F, quote=F)

	#Add gene name using the interval gene correspondence
	interval_2_gene$full_compid <- paste(interval_2_gene$compid,cnv_type_suff,sep="_")
	interval_2_gene_f <- interval_2_gene[,c("V11","full_compid")]
	colnames(interval_2_gene_f)[which(colnames(interval_2_gene_f)=="V11")] <- "gene"
	final_table_fg <- merge(final_table,interval_2_gene_f,by.x="compid",by.y="full_compid",all.x=T)
	write.table(final_table_fg,file=paste(out_f,paste(targetexons_name,out_suff,"final_suspect_table_gene.txt",sep="_"),sep="/"),sep="\t",row.names=F, quote=F)	

#}

	
