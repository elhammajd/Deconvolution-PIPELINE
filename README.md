# Deconvolution-PIPELINE
---

## About PIPELINE

~~~

---
## Directory Layout
![image](image0.png)

We assume the user set the default directory at **beluga** at Compute Canada
~~~
    [PIPELINE]  
~~~
all codes are in the subdirectory directory at **scripts** 
~~~
    [PIPELINE]/scripts  
~~~
all the .sh files that run the R files are in the subdirectory directory at **sh** 
~~~
    [HSPE]/sh  
~~~
all the log files are in the subdirectory directory at **rout** 
~~~
    [PIPELINE]/rout  
~~~
all the final results/intermedia results are in the subdirectory directory at **result** 
~~~
    [PIPELINE]/outputs  
~~~
all the **datasets** from Schelker are stored at the directory bellow, which are accessible to all group members
~~~
    [PIPELINE]/data  
~~~

<details><summary>scripts</summary>

    ├── scripts  
    │ 	 ├── hspe.R		    # install libraries
    │ 	 ├── main.R		# main code to deconvolution
    │ 	 ├── combine_Y_refs.R 			        # function on reference dataset
    │ 	 ├── constr_fns.R			# Normalize matrices and returning predictions  
    │ 	 ├── data.R  	        # funstion on dataset
    │ 	 ├── find_markers_fn.R			# function to build markers    
    │ 	 ├── marker_list.R  	        # function to build marker-list	
    │ 	 ├── process_markers.R  	        # function on processing markers	
    │ 	 └── samplex_sample.R			 				
</details>
<details><summary>sh</summary>

    ├── sh
    │ 	 ├── lib.sh		# sh.file to install the libraries	
    │ 	 └── main.sh		# sh.file to get the deconvolution results				
	
</details>
<details><summary>rout</summary>

    ├── log files after submitting jobs
    │ 	 ├── lib.sh		# log file for lib.sh	
    │ 	 └── main.sh		# log file for main.sh	

</details>
<details><summary>outputs (final & intermedia results)</summary>

    ├──  final result 
    │ 	 ├── predictionsresult.csv			# The final prediction results 
    │ 	 └── plot.png 		# estimade and true proportion
</details>
<details><summary>data</summary>

    ├── data  
    │ 	 ├── scRNAseq_7873_schelker.rds	   
    │ 	 ├── scRNAseq_7882_schelker.rds	
    │ 	 ├── scRNAseq_7892_schelker.rds	
    │ 	 └── bulk-schelker.tsv 		

</details>

---
## Notice

As all the processes are conducted using the relative path, it's very important to set up [PIPELINE] and use it correctly. 
[PIPELINE] should be consisted of three parts: part 1 is ```/project/6003851/``` to ensure all the files can run on Compute Canada; part 2 is your ```user name``` at Compute Canada; part 3 is your ```folder's name```. For example, the writer's directory is as follows:

~~~
/project/6003851/elhma/PIPELINE
~~~

If you are not sure about the path of your working folder, try to type in 'pwd' command in linux or 'getwd()' in R language for reference. 

---
## Before you start
1. decide the path of [PIPELINE] to replicate our results;
2. create the subdirectories **scripts**, **sh**, **rout**, and **outputs** at [PIPELINE]；
3. allocate all relevant files into each subdirectory. The **rout**, and **outputs** folders will be empty at the beginning while the **scripts**, **data** and **sh** folders should look like the figure below:

![image](image1.png)

5. In the main directory [PIPELINE], use the following commands to load R/4.0.2 language in Compute Canada (The environment settings in CC change occasionally, make sure to check and use their latest settings):
~~~
module load gcc/9.3.0 r/4.0.2
mkdir -p ~/.local/R/$EBVERSIONR/
export R_LIBS=~/.local/R/$EBVERSIONR/
~~~
4. before we run the .sh files, we use in the following commands in to install some packages needed for the task
~~~
#packages required in R:
install.packages(c('crayon', 'devtools', 'formatR'))
lapply(pkgs, require, character.only = TRUE)
sessionInfo()
~~~

---

## Running files (estimated time per job)---- 7 days with cpu=64G and ntasks=8   was not enough

<details><summary>1. main.R </summary>

- install lib_hspe;

 </details>
 
 ~~~
    sbatch ./sh/lib.sh
~~~




