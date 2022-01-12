# Compare-TPN-Mandena1
Public code for the manuscript:
Kauffman Kayla, Werner Courtney S., Titcomb Georgia, Pender Michelle, Rabezara Jean Yves, Herrera James P., Shapiro Julie Teresa, Solis Alma, Soarimalala Voahangy, Tortosa Pablo, Kramer Randall, Moody James, Mucha Peter J. and Nunn Charles. 2022 Comparing transmission potential networks based on social network surveys, close contacts and environmental overlap in rural Madagascar. J. R. Soc. Interface. 19:20210690. 20210690. <http://doi.org/10.1098/rsif.2021.0690>

### GPS_data_preparation.R  
(which sources functions in **clean_GPS_functions.R**) goes through the processes:  

- to clean the GPS data  
- select the study area  
- remove data from days when participants likely did not wear a GPS  
- calculate utilization distributions  
- calculate volume intersections  
- calculate proximities  
- check movement fidelity  

### close_contact_predictors.R   
sets up the dataframes for using the in close contact GLM and matricies for the close contact ERGM  

### naming_networks.R  

- builds the directed and undirected networks based off the social network survey data  
- For each network calculate network wide stats and centrality  

### close_contact_network.R  

- builds the ERGM to find the probability of close contact edges  
- builds the GLM to predict the weight of those edges (# contacts / # possible contacts)  
- creates simulated contact networks to calculate network wide stats and centrality on  

### environmental_overlap_networks.R  

- uses classified imagery and individuals' utilization distributions to build the bipartite environmental overlap networks  
- reprojects those networks to a unipartite projection
- calculates network wide stats and centrality

## network_comparisions.R  

All of the network to network comparisions:  

- plots of networks  
- find reciprocated edges on the naming networks  
- edge weight and reciprocated naming comparisons  
    - violin plot comparing full naming and edge weight distributions  
- outlier investigations  
- correlation in centrality  
    - includes "superspreaders" on each network  
