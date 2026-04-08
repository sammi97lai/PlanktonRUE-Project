library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
setwd("~/GitHub/PlanktonRUE-Project/Results/SEM_Figures")

#zooplankton diversity plot ........................................
zoo_div <- grViz("
digraph Phytoplankton_Diversity {

  graph [layout = dot, rankdir =LR, ratio = compress, nodesep = 0.2, ranksep=2]
  node [style=filled, fontname=Helvetica, fontsize=25]
  edge [fontname=Helvetica, fontsize=20, penwidth=3]
  
  Temperature [shape=ellipse, fillcolor=skyblue]
  Salinity [shape=ellipse, fillcolor=skyblue]
  pH [shape=ellipse, fillcolor=skyblue]
  Max_Depth [shape=ellipse, fillcolor=skyblue]
  DOY [shape=ellipse, fillcolor=skyblue]
  Daylight_Hours [shape=ellipse, fillcolor=skyblue]

  Cyanobacteria_Dominance [shape=box, fillcolor=green]
  Phytoplankton_Genus_Richness [shape=box, fillcolor=green]
  Phytoplankton_FEve_selected [shape=box, fillcolor=green]
  Phytoplankton_FDis_selected [shape=box, fillcolor=green]
  
  Zooplankton_Species_Richness [shape=box, fillcolor=magenta]
  Zooplankton_FEve [shape=box, fillcolor=magenta]
  Zooplankton_FDis [shape=box, fillcolor=magenta]

  Salinity -> Zooplankton_Species_Richness [label='β=0.59(0.09,1.09)', color='forestgreen']
  
  Temperature -> Zooplankton_FEve [label='hump-shaped\nβ²=-0.02(-0.03,-0.00)']
  pH -> Zooplankton_FEve [label='β=0.02(0.01,0.04)', color='forestgreen']
  Max_Depth -> Zooplankton_FEve [label='U-shaped\nβ=-0.06(-0.08,-0.04)\nβ²=0.01(0.00,0.02)']
  Phytoplankton_Genus_Richness -> Zooplankton_FEve [label='β=-0.01(-0.01,-0.00)', color='red']
  Phytoplankton_FEve_selected -> Zooplankton_FEve [label='β=-0.12(-0.23,-0.01)', color='red']
  
  Max_Depth -> Zooplankton_FDis [label='U-shaped\nβ=-0.09(-0.10,-0.08)\nβ²=0.02(0.01,0.02)']

}
")
zoo_div

svg <- export_svg(zoo_div)
rsvg_png(charToRaw(svg), "zoo_div_260317.png")

#Zooplankton RUE..............................................................
Zoo_RUE <- grViz("
digraph Zooplankton_RUE {

  graph [layout = dot, rankdir =LR, ratio = compress, nodesep = 0.5]
  node [style=filled, fontname=Helvetica, fontsize=40]
  edge [fontname=Helvetica, fontsize=40, penwidth=3]
  
  Temperature [shape=ellipse, fillcolor=skyblue]
  Salinity [shape=ellipse, fillcolor=skyblue]
  pH [shape=ellipse, fillcolor=skyblue]
  Max_Depth [shape=ellipse, fillcolor=skyblue]
  DOY [shape=ellipse, fillcolor=skyblue]
  Daylight_Hours [shape=ellipse, fillcolor=skyblue]

  Cyanobacteria_Dominance [shape=box, fillcolor=green]
  Phytoplankton_Genus_Richness [shape=box, fillcolor=green]
  Phytoplankton_FEve_selected [shape=box, fillcolor=green]
  Phytoplankton_FDis_selected [shape=box, fillcolor=green]
  
  Zooplankton_Species_Richness [shape=box, fillcolor=magenta]
  Zooplankton_FEve [shape=box, fillcolor=magenta]
  Zooplankton_FDis [shape=box, fillcolor=magenta]

  Zooplankton_RUE [shape=box, fillcolor=orange]
  
  Temperature_and_phyto_genus_richness [shape=box, style='rounded,filled', fillcolor=lightgray]
  Salinity_and_phyto_genus_richness [shape=box, style='rounded,filled', fillcolor=lightgray]
  Temperature_and_Zooplankton_FEve [shape=box, style='rounded,filled', fillcolor=lightgray]
  
# Environment
{ rank = same;
  Temperature; Salinity; pH; Max_Depth; DOY; Daylight_Hours;
}

# Phytoplankton
{ rank = same;
  Cyanobacteria_Dominance;
  Phytoplankton_Genus_Richness;
  Phytoplankton_FEve_selected;
  Phytoplankton_FDis_selected;
  Temperature_and_phyto_genus_richness;
  Salinity_and_phyto_genus_richness;
}

# Zooplankton (NEW separate layer)
{ rank = same;
  Zooplankton_Species_Richness;
  Zooplankton_FEve;
  Zooplankton_FDis;
  Temperature_and_Zooplankton_FEve;
}

# Response
{ rank = same;
  Zooplankton_RUE;
}

# Environment → Phytoplankton
Temperature -> Phytoplankton_Genus_Richness [style=invis]

# Phytoplankton → Zooplankton 
Phytoplankton_Genus_Richness -> Zooplankton_Species_Richness [style=invis]

# Interactions → RUE
Zooplankton_Species_Richness -> Zooplankton_RUE [style=invis]
  
  Temperature -> Zooplankton_RUE [label='β=0.70(0.02,1.39)', color='forestgreen']
  pH -> Zooplankton_RUE [label='β=-0.27(-0.38,-0.16)', color='red']
  Max_Depth -> Zooplankton_RUE [label='β=-0.85(-0.99,-0.72)', color='red']

  Phytoplankton_Genus_Richness -> Zooplankton_RUE [label='β=-0.05(-0.07,-0.03)', color='red']
  Phytoplankton_FEve_selected -> Zooplankton_RUE [label='β=0.20(-0.42,0.83)', color='#90EE90']
  Phytoplankton_FDis_selected -> Zooplankton_RUE [label='β=0.53(-0.64,1.67)', color='#90EE90']
  Cyanobacteria_Dominance -> Zooplankton_RUE [label='β=-1.09(-1.57,-0.61)', color='red']
  Temperature_and_phyto_genus_richness -> Zooplankton_RUE [label='β=-0.03(-0.04,-0.01)', color='red']
  Salinity_and_phyto_genus_richness -> Zooplankton_RUE [label='β=-0.02(-0.05,-0.00)', color='red']
  
  Zooplankton_Species_Richness -> Zooplankton_RUE [label='β=0.03(-0.02,0.07)', color='#90EE90']
  Zooplankton_FEve -> Zooplankton_RUE [label='β=-0.16(-0.56,0.23)', color='#F08080']
  Zooplankton_FDis -> Zooplankton_RUE [label='β=-1.66(-2.37,-0.95)', color='red']
  Temperature_and_Zooplankton_FEve -> Zooplankton_RUE [label='β=0.41(0.02,0.79)', color='forestgreen']
}
")
Zoo_RUE

svg <- export_svg(Zoo_RUE)
rsvg_png(charToRaw(svg), "Zoo_RUE_260317.png")

#all plot......................................................................
Zoo <- grViz("
digraph Zooplankton {

  graph [layout = dot, rankdir =LR, ratio = compress, nodesep = 0.5]
  node [style=filled, fontname=Helvetica, fontsize=40]
  edge [fontname=Helvetica, fontsize=35, penwidth=3]
  
  Temperature [shape=ellipse, fillcolor=skyblue]
  Salinity [shape=ellipse, fillcolor=skyblue]
  pH [shape=ellipse, fillcolor=skyblue]
  Max_Depth [shape=ellipse, fillcolor=skyblue]
  DOY [shape=ellipse, fillcolor=skyblue]
  Daylight_Hours [shape=ellipse, fillcolor=skyblue]

  Cyanobacteria_Dominance [shape=box, fillcolor=green]
  Phytoplankton_Genus_Richness [shape=box, fillcolor=green]
  Phytoplankton_FEve_selected [shape=box, fillcolor=green]
  Phytoplankton_FDis_selected [shape=box, fillcolor=green]
  
  Zooplankton_Species_Richness [shape=box, fillcolor=magenta]
  Zooplankton_FEve [shape=box, fillcolor=magenta]
  Zooplankton_FDis [shape=box, fillcolor=magenta]

  Zooplankton_RUE [shape=box, fillcolor=orange]
  
  Temperature_and_phyto_genus_richness [shape=box, style='rounded,filled', fillcolor=lightgray]
  Salinity_and_phyto_genus_richness [shape=box, style='rounded,filled', fillcolor=lightgray]
  Temperature_and_Zooplankton_FEve [shape=box, style='rounded,filled', fillcolor=lightgray]
  
# Environment
{ rank = same;
  Temperature; Salinity; pH; Max_Depth; DOY; Daylight_Hours;
}

# Phytoplankton
{ rank = same;
  Cyanobacteria_Dominance;
  Phytoplankton_Genus_Richness;
  Phytoplankton_FEve_selected;
  Phytoplankton_FDis_selected;
  Temperature_and_phyto_genus_richness;
  Salinity_and_phyto_genus_richness;
}

# Zooplankton (NEW separate layer)
{ rank = same;
  Zooplankton_Species_Richness;
  Zooplankton_FEve;
  Zooplankton_FDis;
  Temperature_and_Zooplankton_FEve;
}

# Response
{ rank = same;
  Zooplankton_RUE;
}

# Environment → Phytoplankton
Temperature -> Phytoplankton_Genus_Richness [style=invis]

# Phytoplankton → Zooplankton 
Phytoplankton_Genus_Richness -> Zooplankton_Species_Richness [style=invis]

# Interactions → RUE
Zooplankton_Species_Richness -> Zooplankton_RUE [style=invis]
  
  Salinity -> Zooplankton_Species_Richness [label='β=0.59(0.09,1.09)', color='forestgreen']
  
  Temperature -> Zooplankton_FEve [label='hump-shaped\nβ²=-0.02(-0.03,-0.00)']
  pH -> Zooplankton_FEve [label='β=0.02(0.01,0.04)', color='forestgreen']
  Max_Depth -> Zooplankton_FEve [label='U-shaped\nβ=-0.06(-0.08,-0.04)\nβ²=0.01(0.00,0.02)']
  Phytoplankton_Genus_Richness -> Zooplankton_FEve [label='β=-0.01(-0.01,-0.00)', color='red']
  Phytoplankton_FEve_selected -> Zooplankton_FEve [label='β=-0.12(-0.23,-0.01)', color='red']
  
  Max_Depth -> Zooplankton_FDis [label='U-shaped\nβ=-0.09(-0.10,-0.08)\nβ²=0.02(0.01,0.02)']
  
  Temperature -> Zooplankton_RUE [label='β=0.70(0.02,1.39)', color='forestgreen']
  pH -> Zooplankton_RUE [label='β=-0.27(-0.38,-0.16)', color='red']
  Max_Depth -> Zooplankton_RUE [label='β=-0.85(-0.99,-0.72)', color='red']

  Phytoplankton_Genus_Richness -> Zooplankton_RUE [label='β=-0.05(-0.07,-0.03)', color='red']
  Phytoplankton_FEve_selected -> Zooplankton_RUE [label='β=0.20(-0.42,0.83)', color='#90EE90']
  Phytoplankton_FDis_selected -> Zooplankton_RUE [label='β=0.53(-0.64,1.67)', color='#90EE90']
  Cyanobacteria_Dominance -> Zooplankton_RUE [label='β=-1.09(-1.57,-0.61)', color='red']
  Temperature_and_phyto_genus_richness -> Zooplankton_RUE [label='β=-0.03(-0.04,-0.01)', color='red']
  Salinity_and_phyto_genus_richness -> Zooplankton_RUE [label='β=-0.02(-0.05,-0.00)', color='red']
  
  Zooplankton_Species_Richness -> Zooplankton_RUE [label='β=0.03(-0.02,0.07)', color='#90EE90']
  Zooplankton_FEve -> Zooplankton_RUE [label='β=-0.16(-0.56,0.23)', color='#F08080']
  Zooplankton_FDis -> Zooplankton_RUE [label='β=-1.66(-2.37,-0.95)', color='red']
  Temperature_and_Zooplankton_FEve -> Zooplankton_RUE [label='β=0.41(0.02,0.79)', color='forestgreen']
}
")
Zoo

svg <- export_svg(Zoo)
rsvg_png(charToRaw(svg), "Zoo_all_260317.png")
