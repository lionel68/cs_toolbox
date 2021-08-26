---
layout: default
---

**What**: In this workshop we will consider the challenges of different types of CS data ranging from structured, semi-structured and opportunistic data. We will especially focus on accounting for biases (such as spatial) and heterogeneity in observational and sampling processes. 

**For whom**: The course is targeted towards ecologists that:

* have experience using R
* have basic knowledge of Generalized Linear Mixed effect Models
* are working (or planning to work) with CS data


**Who**: Diana Bowler, Lionel Hertzog and Swantje Loebel

**When**: Sunday 29.08.2021, Time 9:30 – 17:30 (CET)

**Where**: Online via [Zoom](https://uni-jena-de.zoom.us/j/61640256162), discussions will be done on [Slack](https://join.slack.com/t/slack-y314603/shared_invite/zt-utiuaxxh-KkJWuTrXMd1KrijG6KT2iw)

### (tentative) Program 

* 09:30 – 09:45: Lecture: Welcoming words, schedule for the session, how the workshop will be organized (chat, repository …) [Swantje]

* 09:45 – 10:15: Lecture: Introduction into the large diversity of citizen science data from structured monitoring programs relying on volunteers (breeding bird survey) to opportunistic observations (GBIF …) and the bias associated with it (site selection, reporting bias, detection issues …) and simple trend models. [lecture file](docs/02_CS_Intro.pdf)

* 10:15 – 10:45: Coding: fitting a trend model (time and space) to structured data, [Rscriptfile](docs/C01_trend_structured.R), [data file 1](docs/C01_temporal_data.rds), [data file 2](docs/C01_spatial_data.rds) 

* 11:00 – 11:30: Lecture: introduction to occupancy models, [lecture file](docs/03_intro_occupancy.pdf)

* 11:30 – 12:00: Coding: fitting occupancy models using R (unmarked) and JAGS, [Rscript file](docs/C02_intro_occupancy.R)

* 12:00 – 13:00: Lunch break

* 13:00 – 13:15: exchange time

* 13:15 – 13:45: Lecture: more advanced occupancy models (more than one species, dynamic occupancy models) [Diana]

* 13:45 – 14:15: Coding: fitting more advanced occupancy models using JAGS [Diana]

* 14:15 – 14:45: Lecture: introduction into spatial models (INLA and inlabru), [lecture file](docs/05_spatial_models.html)

* 14:45 – 15:15: Coding: fitting spatial models (inlabru), [Rscript file](docs/C04_spatial_data.R)

* 15:15 – 15:45: Lecture: Integrated modelling approaches (combined what we learn earlier) [Diana]

* 15:45 – 16:15: Coding: fitting integrated models (INLA, Stan?) [Lionel / Diana]

* 16:15 – End: Open discussion, Stats troubleshooting, data health checks

### Requirements

+ Install [`R`](https://cloud.r-project.org/) and [`RStudio`](https://rstudio.com/products/rstudio/download/#download). 

+ unmarked

+ INLA

+ inlabru

+ rjags

+ tidyverse

### To-do list

+ ...
