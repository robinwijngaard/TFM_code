PattRec - CNV detection tool
=============================


DESCRIPTION
------------
PattRec is a bioinformatics tool designed to detect rare copy number variants (CNVs) in targeted Next Generation Sequencing (tg-NGS) data. It is presented as a Java-based GUI, with its CNV detection algorithm implemented in R.
This tool was designed for use with target gene panels, sequenced in Illumina platforms.




SYSTEM REQUIREMENTS
------------
- **Ubuntu 16.04 LTS** operating system or greater.
- Java 8.
- MySQL or MariaDB.
- 8GB RAM or greater is recommended.



INSTALLATION
------------
This application must be installed from command line.
To install it, execute the following commands:

1. If you don't have Java 8 installed, first you need to add Java PPA repository in your system to install Oracle Java 8 using the following commands:
```
sudo add-apt-repository ppa:webupd8team/java
sudo apt-get update
```

2. Move to the directory where you downloaded pattrec.deb and type the following commands:
```
sudo dpkg -i pattrec.deb
```
This command will unpack the package. Probably you are going to get an error of missing dependencies, then execute the following:
```
sudo apt-get install -f
```
This second command will fix the problem and install missing dependencies and follow up the installation of pattrec.

If your system didn't have yet Java 8 installed, a prompt will be shown during the installation to accept Oracle license.
If you haven't any database installed in your system, you will be asked for a password for your root user. You could use this later in PattRec for the database configuration process.



EXECUTION
------------
To execute PattRec, you can search it in your system applications, or you can type on a terminal:
```
pattrec
```
The first time the application is executed, you will be asked for dataset configuration. You can introduce the root user you have configured during installation, or any other user you have to access MySQL database.
The application will create a folder named 'PattRec' in your home directory. The results of the executions will be placed there.


TROUBLESHOOTING
------------
###### MySQL root password missing.
If you did not have a previous installation of MySQL or you do not have a user to manage the database in PattRec, you can create one in the following way: 
1. Enter the database with the following bash command:
```
sudo mysql
```
2. Create a new user and grant privileges with the following mysql commands:
```
CREATE USER 'user_name'@'localhost' IDENTIFIED BY 'password';
GRANT ALL_PRIVILEGES ON *.* TO 'user_name'@'localhost' WITH GRANT OPTION;
FLUSH PRIVILEGES;
```
Now you can use these credentials to configure the database in PattRec.


USER MANUAL
------------
1. Configure database: click "Database configuration" and enter user and password for MySQL or MariaDB database.

2. Select input files:
	- BAM test: BAM file of the test sample (patient).
	- BAM control: one or more BAM files to compare with the BAM test.
	
		*NOTE: It is highly recommended that all the BAM files were sequenced on the same run.*
	- BED file: file containing the target regions (sequenced regions) in BED format.
	- FASTA file: genome reference file (it must be the same used for the alignment of the BAM files).
		*NOTE: GRCh37-hg19 can be downloaded from: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/*
		*GRCh38-hg38 can be downloaded from: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/*

3. Click "RUN".

4. Configure parameters for the algorithm (if needed).

5. Click "Proceed".

6. You will be asked to save results in database. If you agree, that results will be used in the next execution of the algorithm.

Results will be placed in the user home, in the folder "PattRec".



