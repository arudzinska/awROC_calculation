# awROC_calculation

Module to generate Receiver Operating Characteristic (**ROC**) curves, calculate **ROCE** (ROC Enrichments) and **AUC** (Area Under the Curve) - the standard ones, as well as their average weighted modification (**awROC** curve, **awROCE**, **awAUC**) which includes active compounds' clustering information. These metrics are used in virtual screening tools benchmarking tests to assess the performance of the software.

The module reads the output virtual screening ranking from PharmScreen, Pharmacelera's tool for ligand-based virtual screening (see [the webpage of Pharmacelera](https://www.pharmacelera.com/)). The output of the calculation is a CSV file with the enrichments (at 0.5%, 1%, 2% and 5% of false positives fraction retrieved), AUC and a PNG file with the ROC curve.

### ROC
ROC curve renders the ability of the tool to distinguish between two populations: true active compounds and decoys - inactive molecules. X and Y values of the ROC curve at the given point are calculated as follows:

![equation-awroc](https://lh6.googleusercontent.com/3Kjk2NAp9RFM6QeLDj8vIcV2rGwIjJzC3yJcctEFWlrE2tjjAigniNAKcaFuFbUg8SskbC41NfEwVBQafKH-=w1366-h635)

where: *X%* is the fraction of the decoys retrieved at the chosen position of the virtual screening ranking.

When dividing the Y point value by the X point value one obtains the ROC Enrichment at the given retrieved decoys fraction. AUC is the area under the whole ROC curve.

### awROC
The average weighted modification inlcudes information about active compounds' clustering to evaluate the tool's ability to retrieve new scaffolds. The modified equation for awROC curve points and awROC Enrichments is as follows:

![equation-awroce](https://lh3.googleusercontent.com/BDZF94Qp6gvr9exJZ6_hnjuyFakQpBDsk3YQaQSmgyKJfMQ8d5ceUHd34daasVrTQmSWOA305UTFy9rbZAne=w1366-h635)

where: *w<sub>ij</sub>* = 1/*N<sub>j</sub>* and is the weight of the *i*th structure from the *j*th cluster. *N<sub>j</sub>* is the number of structures in given cluster. *&alpha;<sup>X%</sup><sub>ij</sub>* is 1 or 0 depending on whether the *i*th structure of the *j*th cluster already (respectively) appeared or not in the chosen fraction of the dataset.

Similarly to the standard ROC curve, the awROC Enrichment can be calculated by dividing the Y point value by the X point value of the curve and awAUC is simply the area under the obtained curve.
