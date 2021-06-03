# Draw cell lineage tree with gene expression
This program can draw cell lineage tree with gene expression into svg file from one or two CD files.

<img src="https://github.com/xthuang226/Draw_cell_lineage_tree/blob/main/results/CD130826PHA4p2_CD140102NHR25p2_blot_129t_349c_P0.svg">

## Here are all options can be used:  

### [Input and Output]  
--ref  "Specify the reference CD file."  
--align  "Specify the CD file to be alinged."  
-o, --output  default='results'  "Specify the output file folder"

### [Control the data to be draw on the tree]  
-c, --column  default='blot'  "Select a specific column to show on the tree."  
-r, --root  default='P0'  "Specify the root of the tree."  
--endtp  "Set an end time point of the tree."  
--cellstage  "Set a cellstage of the tree."  

### [Control other drawings except the tree]  
-nl, --nolabel  default=True  "Do not draw cell label."  
-na, --noaxis  default=True  "Do not draw the time axis."  
-nb, --nobrand  default=True  "Do not draw the brand."  
-nt, --notitle  default=True  "Do not draw the title."  

### [Control the drawing line]  
--linewidth  default=5  "Set the linewidth."  
--lineinter  default=15  "Set the lineinter."  

### [Control the color]  
--clvl  default=10  "The number of color levels to be draw."  
--clvlmax  default=0  "The number of color levels in the palette."  
--clvllow  default=1, help="Low brightness part startting position."  
--clvllowplus  default=1  "The number of levels in the low brightness part."  
--clvlhigh  default=0  "High brightness part startting position."  
--revc  default=0  "Reverse the reference and the align color."  
--cl1  default=255  "The first color max brightness."  
--cl2  default=255  "The second color max brightness."  
--seaborn  default=0  "Use seaborn color palette."  


## Here are some examples of running the program:  

### [Draw tree with one CD file]  
python Draw_cell_lineage_tree.py --ref data/CD130826PHA4p2.csv --cellstage 350  
python Draw_cell_lineage_tree.py --ref data/CD130826PHA4p2.csv --endtp 120  

### [Draw tree with two CD files]  
python Draw_cell_lineage_tree.py --ref data/CD130826PHA4p2.csv --align data/CD140102NHR25p2.csv --cellstage 350  

### [Others]  
python Draw_cell_lineage_tree.py --ref data/CD130826PHA4p2.csv --align data/CD140102NHR25p2.csv --cellstage 350 -r ABarp  
python Draw_cell_lineage_tree.py --ref data/CD130826PHA4p2.csv --align data/CD140102NHR25p2.csv --cellstage 350 -c none  

