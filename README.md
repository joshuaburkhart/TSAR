#TSAR: Time-Series Analysis tool for Respiratory Viral DREAM Challenge (syn7202976)
Joshua Burkhart, M.Sc.
Department of Medical Informatics and Clinical Epidemiology
Oregon Health & Science University
Portland, Oregon 97239

##This project uses supervised microarray normalization and support vector machine regression to predict log symptom scores from probe intensity for the Respiratory Viral DREAM Challenge[0].

##Introduction
The provided RMA-normalized microarray data was found to suffer from study-specific batch effects which were corrected for using snm[1]. The probes whose intensity varied most among subjects' snm-normalized microarray data between earliest and nearest to 24h assays were selected as features and used to train and tune e1071 support vector machines[2].

&nbsp;

Supervised normalization was used after finding study-specific batch effects among the provided microarray data. Additionally, in order to attempt to find a sex-inspecific signal, supervised normalization was used to adjust for putative gender effects.

&nbsp;

Due to the robustness and tuning capability of the e1071 R library, support vector machine regression was chosen for modeling.

&nbsp;

Features were selected following supervised normalization of microarray data. Because the clinical training data indicated subjects tended to report more symptoms at 0h, it was thought a nocebo-effect could irregularly alter later expression data. Therefore, only the earliest microarray data was used for each subject in generating the 0h predictions. The symptom score training data showed subjects reported the most symptoms following 24h. Assuming true gene expression signal would be well correlated with reported symptoms, it was thought the microarray data closest to (but not beyond) 24h would contain the highest signal to noise ratio relative to the other timepoints.

##Methods
###Common between 0h and 24h
The provided RMA expression and clinical training data was filtered to remove sham subjects and unused timepoints. Only each subject's earliest and latest timepoints (<= 24h) were retained. Supervised normalization was applied to the remaining expression data using timepoint (earliest or latest) as the biological variables and both study and gender as adjustment variables. The 2000 most variable non-control probes following supervised normalization were retained as features. 

###0h
Under the assumption that non-logarithmic scores may be easier to learn, each training subject's log symptom score was transformed (10^LOGSYMPTSCORE_SC3^) to yield a symptom score. Support vector machines were tuned using k=10 fold cross-validation on each subject's earliest snm-normalized microarray data varying epsilon (from 0 to 1 by increments of 0.1) and cost (from 2^2^ to 2^9^ by increments of powers of two) to predict the symptom score. The best model found using k=10 fold cross-validation was used to perform leave-one-out cross-validation (LOOCV) on each subject's earliest snm-normalized microarray data (predicting symptom scores). The LOOCV predictions were then log10-transformed to yield a predicted LOGSYMPTSCORE_SC3 value.

&nbsp;

The provided RMA-normalized and clinical test data (Phase 1, Phase 2, and Phase 3) was filtered to remove sham subjects (0) and unused timepoints. Only each subject's earliest timepoint was retained. The best model found during training was used to predict the symptom score for each subject, which was then log10-transformed to yield a predicted LOGSYMPTSCORE_SC3 value.

###24h
Similarly to the 0h prediction solution outlined above, each training subject's log symptom score was transformed (10^LOGSYMPTSCORE_SC3^) to yield a symptom score. Support vector machines were tuned using k=10 fold cross-validation on each subject's *latest* snm-normalized microarray data varying epsilon (from 0 to 1 by increments of 0.1) and cost (from 2^2^ to 2^9^ by increments of powers of two) to predict the symptom score. The best model found using k=10 fold cross-validation was used to perform leave-one-out cross-validation (LOOCV) on each subject's *latest* snm-normalized microarray data (predicting symptom scores). The LOOCV predictions were then log10-transformed to yield a predicted LOGSYMPTSCORE_SC3 value.

&nbsp;

Also similarly to the 0h prediction solution, the provided RMA-normalized and clinical test data (Phase 1, Phase 2, and Phase 3) was filtered to remove sham subjects (0) and unused timepoints. Only each subject's *latest* timepoint was retained. The best model found during training was used to predict the symptom score for each subject, which was then log10-transformed to yield a predicted LOGSYMPTSCORE_SC3 value.

##Conclusion/Discussion
Though supervised microarray normalization appeared to be an effective remedy, the study-specific batch effects present in this challenge may limit the power to discover an expression signature well correlated with symptom scores. Additionally, several subjects reported symptoms at or shortly after 0h, indicating a nocebo-effect was present. This could perhaps be ameliorated some by more objective symptom quantification. 

## Submission Requirements
Expression Time Range | Predictors | Code | Leave-one-out CVs | Official Synapse Submission File | Official Synapse Submission ID
-|-|-|-|-|-
`Up to Hour 0`| syn7207901 | syn7207894 | syn7207900 | syn7207899 | 7207905
`Up to Hour 24`| syn7207898 | syn7207896 | syn7207897 | syn7207895 | 7207906

##References
[0] Respiratory Viral DREAM Challenge (syn5647810)
[1] Mecham BH, Nelson PS and Storey JD (2010). “Supervised normalization of microarrays.” Bioinformatics, 26, pp. 1308-1315. 
[2] Meyer, David, et al. "Package ‘e1071’." CRAN R Project (2015).
[3] TSAR: Time-Series Analysis tool for Respiratory Viral DREAM Challenge (syn7202976)
