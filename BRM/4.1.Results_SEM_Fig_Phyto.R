library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
setwd("~/GitHub/PlanktonRUE-Project/Figures/SEM_Figures")

#phytoplankton diversity plot ........................................
phyto_div <- grViz("
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
  Genus_Richness [shape=box, fillcolor=green]
  FEve [shape=box, fillcolor=green]
  FDis [shape=box, fillcolor=green]

  Salinity -> Cyanobacteria_Dominance [label='β=-2.88(-3.93,-1.81)\nβ²=-0.85(-1.24,-0.46)', color='forestgreen']
  pH -> Cyanobacteria_Dominance [label='U-shaped\nβ²=0.07(0.03,0.11)']
  Max_Depth -> Cyanobacteria_Dominance [label='U-shaped\nβ²=0.04(0.00,0.07)']
  
  Temperature -> Genus_Richness [label='β=-2.15(-2.94,-1.37)',color='red']
  Salinity -> Genus_Richness [label='U-shaped\nβ²=1.09(0.71,1.47)']
  pH -> Genus_Richness [label='β=0.84(0.51,1.17)\nβ²=-0.15(-0.30,-0.00)',color='forestgreen']
  Max_Depth -> Genus_Richness [label='β=0.70(0.16,1.26)',color='forestgreen']
  DOY -> Genus_Richness [label='hump-shaped\nβ²=-0.84(-1.44,-0.22)']
  
  Salinity -> FEve [label='U-shaped\nβ²=0.02(0.01,0.03)']
  
  Temperature -> FDis [label='U-shaped\nβ²=0.00(0.00,0.01)']
  pH -> FDis [label='U-shaped\nβ²=-0.01(-0.01,-0.00)']
  Max_Depth -> FDis [label='β=0.01(0.00,0.01)',color='forestgreen']
  DOY -> FDis [label='hump-shaped\nβ²=-0.04(-0.07,-0.00)']
}
")
phyto_div

svg <- export_svg(phyto_div)
rsvg_png(charToRaw(svg), "phyto_div_260317.png")

#RUE per P plot............................................................
RUE_P <- grViz("
digraph Phytoplankton_RUE {

  graph [layout = dot, rankdir =LR, ratio = compress, nodesep = 0.2]
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
  Phytoplankton_FEve [shape=box, fillcolor=green]
  Phytoplankton_FDis [shape=box, fillcolor=green]
  
  Zooplankton_Species_Richness [shape=box, fillcolor=magenta]
  Zooplankton_FEve [shape=box, fillcolor=magenta]
  Zooplankton_FDis [shape=box, fillcolor=magenta]

  RUE_per_P [shape=box, fillcolor=orange]
  
  Salinity_and_Phytoplankton_FDis [shape=box, style='rounded,filled', fillcolor=lightgray]
  Salinity_and_Cyanobacteria_Dominance [shape=box, style='rounded,filled', fillcolor=lightgray]
  Temperature_and_Zooplankton_FEve [shape=box, style='rounded,filled', fillcolor=lightgray]
  Temperature_and_Zooplankton_FDis [shape=box, style='rounded,filled', fillcolor=lightgray]
  
# Environment
{ rank = same;
  Temperature; Salinity; pH; Max_Depth; DOY; Daylight_Hours;
}

# Phytoplankton
{ rank = same;
  Cyanobacteria_Dominance;
  Phytoplankton_Genus_Richness;
  Phytoplankton_FEve;
  Phytoplankton_FDis;
  Salinity_and_Phytoplankton_FDis;
  Salinity_and_Cyanobacteria_Dominance;
}

# Zooplankton (NEW separate layer)
{ rank = same;
  Zooplankton_Species_Richness;
  Zooplankton_FEve;
  Zooplankton_FDis;
  Temperature_and_Zooplankton_FEve;
  Temperature_and_Zooplankton_FDis;
}

# Response
{ rank = same;
  RUE_per_P;
}

# Environment → Phytoplankton
Temperature -> Phytoplankton_Genus_Richness [style=invis]

# Phytoplankton → Zooplankton 
Phytoplankton_Genus_Richness -> Zooplankton_Species_Richness [style=invis]

# Interactions → RUE
Zooplankton_Species_Richness -> RUE_per_P [style=invis]
  
  Temperature -> RUE_per_P [label='U-shaped\nβ²=0.04(0.01,0.08)']
  Salinity -> RUE_per_P [label='β=-0.84(-1.25,-0.45)',color='red']
  pH -> RUE_per_P [label='β=0.13(0.08,0.17)',color='forestgreen']
  DOY -> RUE_per_P [label='hump-shaped\nβ=-2.88(-3.93,-1.81)\nβ²=-0.85(-1.24,-0.46)']
  Daylight_Hours -> RUE_per_P [label='β=-1.90(-2.61,-1.18)\nβ²=0.22(0.11,0.32)',color='red']
  
  Phytoplankton_Genus_Richness -> RUE_per_P [label='β=0.04(0.03,0.05)',color='forestgreen']
  Phytoplankton_FEve -> RUE_per_P [label='β=-0.69(-0.97,-0.41)',color='red']
  Phytoplankton_FDis -> RUE_per_P [label='β=-2.31(-3.06,-1.55)',color='red']
  Cyanobacteria_Dominance -> RUE_per_P [label='β=1.10(0.88,1.31)',color='forestgreen']
  
  Zooplankton_Species_Richness -> RUE_per_P [label='β=0.01(-0.01,0.03)',color='#90EE90']
  Zooplankton_FEve -> RUE_per_P [label='β=-0.05(-0.21,0.11)',color='lightred']
  Zooplankton_FDis -> RUE_per_P [label='β=0.42(0.14,0.70)',color='forestgreen']
  
  Salinity_and_Phytoplankton_FDis -> RUE_per_P [label='β=4.27(3.46,5.07)',color='forestgreen']
  Salinity_and_Cyanobacteria_Dominance -> RUE_per_P [label='β=-0.34(-0.62,-0.07)',color='red']
  Temperature_and_Zooplankton_FEve -> RUE_per_P [label='β=-0.21(-0.38,-0.04)',color='red']
  Temperature_and_Zooplankton_FDis -> RUE_per_P [label='β=0.41(0.08,0.73)',color='forestgreen']
  
}
")
RUE_P

svg <- export_svg(RUE_P)
rsvg_png(charToRaw(svg), "RUE_P_260317.png")

#RUE per N plot.................................................................
RUE_N <- grViz("
digraph Phytoplankton_RUE_per_N {

  graph [layout = dot, rankdir =LR, ratio = compress, nodesep = 0.2]
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
  Phytoplankton_FEve [shape=box, fillcolor=green]
  Phytoplankton_FDis [shape=box, fillcolor=green]
  
  Zooplankton_Species_Richness [shape=box, fillcolor=magenta]
  Zooplankton_FEve [shape=box, fillcolor=magenta]
  Zooplankton_FDis [shape=box, fillcolor=magenta]

  RUE_per_N [shape=box, fillcolor=orange]
  
  Salinity_and_Phytoplankton_FDis [shape=box, style='rounded,filled', fillcolor=lightgray]
  Temperature_and_Phytoplankton_FEve [shape=box, style='rounded,filled', fillcolor=lightgray]
  Salinity_and_Cyanobacteria_Dominance [shape=box, style='rounded,filled', fillcolor=lightgray]
  
# Environment
{ rank = same;
  Temperature; Salinity; pH; Max_Depth; DOY; Daylight_Hours;
}

# Phytoplankton
{ rank = same;
  Cyanobacteria_Dominance;
  Phytoplankton_Genus_Richness;
  Phytoplankton_FEve;
  Phytoplankton_FDis;
  Salinity_and_Phytoplankton_FDis;
  Salinity_and_Cyanobacteria_Dominance;
  Temperature_and_Phytoplankton_FEve;
}

# Zooplankton (NEW separate layer)
{ rank = same;
  Zooplankton_Species_Richness;
  Zooplankton_FEve;
  Zooplankton_FDis;
}

# Response
{ rank = same;
  RUE_per_N;
}

# Environment → Phytoplankton
Temperature -> Phytoplankton_Genus_Richness [style=invis]

# Phytoplankton → Zooplankton 
Phytoplankton_Genus_Richness -> Zooplankton_Species_Richness [style=invis]

# Interactions → RUE
Zooplankton_Species_Richness -> RUE_per_N [style=invis]
  
  Temperature -> RUE_per_N [label='U-shaped\nβ²=0.08(0.04,0.11)']
  Salinity -> RUE_per_N [label='β=-0.56(-0.91,-0.21)',color='red']
  pH -> RUE_per_N [label='β=0.15(0.10,0.19)',color='forestgreen']
  Max_Depth -> RUE_per_N [label='β=-0.07(-0.14,-0.01)',color='red']
  DOY -> RUE_per_N [label='hump-shaped\nβ=-2.23(-3.17,-1.24)\nβ²=-0.68(-1.02,-0.31)']
  Daylight_Hours -> RUE_per_N [label='β=-1.46(-2.09,-0.79)\nβ²=0.16(0.06,0.26)',color='red']
  
  Phytoplankton_Genus_Richness -> RUE_per_N [label='β=0.03(0.03,0.04)',color='forestgreen']
  Phytoplankton_FEve -> RUE_per_N [label='β=-0.84(-1.09,-0.58)',color='red']
  Phytoplankton_FDis -> RUE_per_N [label='β=-2.49(-3.17,-1.82)',color='red']
  Cyanobacteria_Dominance -> RUE_per_N [label='β=1.10(0.90,1.29)',color='forestgreen']
  
  Zooplankton_Species_Richness -> RUE_per_N [label='β=0.01(-0.01,0.02)',color='#90EE90']
  Zooplankton_FEve -> RUE_per_N [label='β=-0.08(-0.23,0.07)',color='#F08080']
  Zooplankton_FDis -> RUE_per_N [label='β=0.17(-0.08,0.42)',color='#90EE90']
  
  Temperature_and_Phytoplankton_FEve -> RUE_per_N [label='β=0.28(0.02,0.54)',color='forestgreen']
  Salinity_and_Phytoplankton_FDis -> RUE_per_N [label='β=4.25(3.50,4.99)',color='forestgreen']
  Salinity_and_Cyanobacteria_Dominance -> RUE_per_N [label='β=-0.30(-0.54,-0.05)',color='red']
}
")
RUE_N

svg <- export_svg(RUE_N)
rsvg_png(charToRaw(svg), "RUE_N_260317.png")

#all plot.........................................................................
grViz("
digraph Phytoplankton_RUE {

  graph [layout = dot, rankdir = LR, splines = true]
  node [style=filled]
  
  Temperature [shape=ellipse, fillcolor=skyblue]
  Salinity [shape=ellipse, fillcolor=skyblue]
  pH [shape=ellipse, fillcolor=skyblue]
  Max_Depth [shape=ellipse, fillcolor=skyblue]
  DOY [shape=ellipse, fillcolor=skyblue]
  Daylight_Hours [shape=ellipse, fillcolor=skyblue]

  Cyanobacteria_Dominance [shape=box, fillcolor=green]
  Phytoplankton_Genus_Richness [shape=box, fillcolor=green]
  Phytoplankton_FEve [shape=box, fillcolor=green]
  Phytoplankton_FDis [shape=box, fillcolor=green]
  
  Zooplankton_Species_Richness [shape=box, fillcolor=magenta]
  Zooplankton_FEve [shape=box, fillcolor=magenta]
  Zooplankton_FDis [shape=box, fillcolor=magenta]

  RUE_per_P [shape=box, fillcolor=orange]
  RUE_per_N [shape=box, fillcolor=orange]
  
  Salinity_and_Phytoplankton_FDis [shape=box, fillcolor=white]
  Temperature_and_Phytoplankton_FEve [shape=box, fillcolor=white]
  Salinity_and_Cyanobacteria_Dominance [shape=box, fillcolor=white]
  Temperature_and_Zooplankton_FEve [shape=box, fillcolor=white]
  Temperature_and_Zooplankton_FDis [shape=box, fillcolor=white]
  
  # Create a “dummy” left → middle → right chain
Temperature -> Phytoplankton_Genus_Richness [style=invis]
Salinity -> Phytoplankton_Genus_Richness [style=invis]
pH -> Phytoplankton_Genus_Richness [style=invis]
Max_Depth -> Phytoplankton_Genus_Richness [style=invis]
DOY -> Phytoplankton_Genus_Richness [style=invis]
Daylight_Hours -> Phytoplankton_Genus_Richness [style=invis]

# Optional: middle → RUE invisible edges to preserve ordering
Phytoplankton_Genus_Richness -> RUE_per_P [style=invis]
Phytoplankton_FEve -> RUE_per_P [style=invis]
Phytoplankton_FDis -> RUE_per_P [style=invis]
Cyanobacteria_Dominance -> RUE_per_P [style=invis]
Zooplankton_FEve -> RUE_per_P [style=invis]
Zooplankton_FDis -> RUE_per_P [style=invis]
Zooplankton_Species_Richness -> RUE_per_P [style=invis]
  
  Salinity -> Cyanobacteria_Dominance [label='β=-2.88(-3.93,-1.81),β²=-0.85(-1.24,-0.46)', color='forestgreen']
  pH -> Cyanobacteria_Dominance [label='U-shaped, β²=0.07(0.03,0.11)']
  Max_Depth -> Cyanobacteria_Dominance [label='U-shaped, β²=0.04(0.00,0.07)']
  
  Temperature -> Phytoplankton_Genus_Richness [label='β=-2.15(-2.94,-1.37)',color='red']
  Salinity -> Phytoplankton_Genus_Richness [label='U-shaped, β²=1.09(0.71,1.47)']
  pH -> Phytoplankton_Genus_Richness [label='β=0.84(0.51,1.17),β²=-0.15(-0.30,-0.00)',color='forestgreen']
  Max_Depth -> Phytoplankton_Genus_Richness [label='β=0.70(0.16,1.26)',color='forestgreen']
  DOY -> Phytoplankton_Genus_Richness [label='hump-shaped, β²=-0.84(-1.44,-0.22)']
  
  Salinity -> Phytoplankton_FEve [label='U-shaped, β²=0.02(0.01,0.03)']
  
  Temperature -> Phytoplankton_FDis [label='U-shaped, β²=0.00(0.00,0.01)']
  pH -> Phytoplankton_FDis [label='U-shaped, β²=-0.01(-0.01,-0.00)']
  Max_Depth -> Phytoplankton_FDis [label='β=0.01(0.00,0.01)',color='forestgreen']
  DOY -> Phytoplankton_FDis [label='hump-shaped, β²=-0.04(-0.07,-0.00)']
  
  Temperature -> RUE_per_P [label='U-shaped, β²=0.04(0.01,0.08)']
  Salinity -> RUE_per_P [label='β=-0.84(-1.25,-0.45)',color='red']
  pH -> RUE_per_P [label='β=0.13(0.08,0.17)',color='forestgreen']
  DOY -> RUE_per_P [label='hump-shaped, β=-2.88(-3.93,-1.81),β²=-0.85(-1.24,-0.46)']
  Daylight_Hours -> RUE_per_P [label='β=-1.90(-2.61,-1.18),β²=0.22(0.11,0.32)',color='red']
  
  Phytoplankton_Genus_Richness -> RUE_per_P [label='β=0.04(0.03,0.05)',color='forestgreen']
  Phytoplankton_FEve -> RUE_per_P [label='β=-0.69(-0.97,-0.41)',color='red']
  Phytoplankton_FDis -> RUE_per_P [label='β=-2.31(-3.06,-1.55)',color='red']
  Cyanobacteria_Dominance -> RUE_per_P [label='β=1.10(0.88,1.31)',color='forestgreen']
  
  Zooplankton_FDis -> RUE_per_P [label='β=0.42(0.14,0.70)',color='forestgreen']
  
  Salinity_and_Phytoplankton_FDis -> RUE_per_P [label='β=4.27(3.46,5.07)',color='forestgreen']
  Salinity_and_Cyanobacteria_Dominance -> RUE_per_P [label='β=-0.34(-0.62,-0.07)',color='red']
  Temperature_and_Zooplankton_FEve -> RUE_per_P [label='β=-0.21(-0.38,-0.04)',color='red']
  Temperature_and_Zooplankton_FDis -> RUE_per_P [label='β=0.41(0.08,0.73)',color='forestgreen']
  
  Temperature -> RUE_per_N [label='U-shaped, β²=0.08(0.04,0.11)']
  Salinity -> RUE_per_N [label='β=-0.56(-0.91,-0.21)',color='red']
  pH -> RUE_per_N [label='β=0.15(0.10,0.19)',color='forestgreen']
  Max_Depth -> RUE_per_N [label='β=-0.07(-0.14,-0.01)',color='red']
  DOY -> RUE_per_N [label='hump-shaped, β=-2.23(-3.17,-1.24),β²=-0.68(-1.02,-0.31)']
  Daylight_Hours -> RUE_per_N [label='β=-1.46(-2.09,-0.79),β²=0.16(0.06,0.26)',color='red']
  
  Phytoplankton_Genus_Richness -> RUE_per_N [label='β=0.03(0.03,0.04)',color='forestgreen']
  Phytoplankton_FEve -> RUE_per_N [label='β=-0.84(-1.09,-0.58)',color='red']
  Phytoplankton_FDis -> RUE_per_N [label='β=-2.49(-3.17,-1.82)',color='red']
  Cyanobacteria_Dominance -> RUE_per_N [label='β=1.10(0.90,1.29)',color='forestgreen']
  
  Temperature_and_Phytoplankton_FEve -> RUE_per_N [label='β=0.28(0.02,0.54)',color='forestgreen']
  Salinity_and_Phytoplankton_FDis -> RUE_per_N [label='β=4.25(3.50,4.99)',color='forestgreen']
  Salinity_and_Cyanobacteria_Dominance -> RUE_per_N [label='β=-0.30(-0.54,-0.05)',color='red']
}
")
