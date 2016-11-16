# Clean-and-Combine
Clean and Combine is currently a package I am building which makes life a little easier in metabolomics data processing. 
Often a lab will look at the same targeted metabolite profile in samples.  In some mass spec integration and quantification softwares, after integrating the peaks, the targeted metabolite data is output as text file which needs to be copied to excel. 
From this the user must then remove concentrations or values below the limit of detection (LOD) or limit of quantification (LOQ),
then merge all of the targeted metabolites into a single data sheet.  This is not only time consuming but leaves room for human error to
compromise the integrity of the data.

The goal of Clean and Combine is to automate several of these steps we often face in targeted metabolomics to make life easier. 
Currently it just a programmable script, where you can use find and replace to change the name of the metabolite to help streamline 
the process, by highlights from the silenced text (i.e., ######) to the last line which removes all R objects to clean up the environment. 
In my example linoleic acid (LA) is the first metabolite to be imported.  I usually use TargetLynx (Waters software) so that is what my
example is based on.  After importing the data set for LA, we make the rownames the files names which will help us when merging all the
metabolites data together. Next we remove blanks from the data set using grepl, and create objects that will serve as our limit of 
detection (LOD) and limit of quantification (LOQ) flags.  The standards are then split from the data set, where we flag the max 
concentration of the standard, which we can use to evaluate the number of samples that are above our standard curve. This can be
visualized using the plot, where the retention time is displayed on the x-axis and the concentration on the y-axis, which can help
you see whether there are any outliers based on a possible wrong peak integration.  Also, the standards are highlighted using a
different colored point.

Next, the concentration is normalized based on your dilution factors etc.  The data frame is then saved as a unique object that can be
called up later during data processing if further investigation into a metabolite is needed. The analytes are split from the data set,
then in my example they are compiled into a file called Ox.dat.  If the rownames of the new data frame do not match that of the analytes
an error should appear. A set of dataframes are then generated will contain information on the number of samples above your highest 
standard curve concentration, the minimum detectable value which you can later use to impute values (i.e., sqrt(minimum detectable value),
the number of samples below the LOQ, and then these dataframes are saved as objects. All R objects and data frames are then removed, since
all information is saved. 

Now if we were import the data for 13_hode, the same sequence of steps would be applied.  The only difference is the information based on 
the code under 13_hode would be complied with that from LA.  All the R data frames that were saved are loaded back and the information 
from 13_hde is added, so the user can view information on all the metabolites in a unique dataframe (i.e., all oxylipid data [Ox.dat],
the number of samples above the highest standard curve for each oxylipid [Max.val], the minimum detectable value for each oxylipid 
[Min.val], etc)

In closing, the first metabolite that you want to enter your dataframe, you will modify the code for LA.  For all metabolites, you want
to enter your dataframe after your first metabolite, you use the code for 13_hode.  The name of the metabolites in all instances can be
easily changed by using find and replace. I working on this as a package after I finish my PhD, but this has been extremely helpful for
cleaning and combing data from single csv files in targeted metabolomic analyses.
