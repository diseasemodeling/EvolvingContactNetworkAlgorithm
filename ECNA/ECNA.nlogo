extensions[
  nw
  csv
  matrix
]


breed [infecteds infectedPeople]
breed [susceptibles susceptible]
breed [non-agent-susceptibles NAsusceptible]

undirected-link-breed [ECNA-links ECNA-link]

turtles-own [
  degree ;Degree or number of contacts, num-IDU-partner
  clustering ;Clustering coefficient
  number-exposures ;Number exposures at a time step, IDU-acts-monthly-counter
  time-in-simulation ;Time turtles has been in the simulation
  desired-degree
  stage ;will get rid of the need for turtles with different breeds
  time-of-infection
  time-of-diagnosis
  source-of-infection
  ID
]
infecteds-own [
  aware? ;Boolean, Is agents aware of their HIV positive status, aware?
  time-discovered ;How long for agent to become diagnosed, age-Diag
  time-infected ;At what point did the agent become infected, age-at-infection

  degree ;Degree or number of contacts, num-IDU-partner
  clustering ;Clustering coefficient
  number-exposures ;Number exposures at a time step, IDU-acts-monthly-counter
  time-in-simulation ;Time turtles has been in the simulation
  desired-degree
  stage ;will get rid of the need for turtles with different breeds
  time-of-infection
  time-of-diagnosis
  source-of-infection
  ID
  partners-per-month
]

susceptibles-own [
  eligible? ;Boolean, Is this susceptible agent an eligible contact for all newly-infected turtle
  used? ;Boolean, Is this susceptible agent used in ERGM
  contacts? ;Does not influence results, measures the number of agents that could possibly get infected at a time-step
  pot-triad? ;Delete, True if a newly infected agents could create a triad with this susceptible contact

  degree ;Degree or number of contacts, num-IDU-partner
  clustering ;Clustering coefficient
  number-exposures ;Number exposures at a time step, IDU-acts-monthly-counter
  time-in-simulation ;Time turtles has been in the simulation
  desired-degree
  stage ;will get rid of the need for turtles with different breeds
  time-of-infection
  time-of-diagnosis
  source-of-infection
  ID
  partners-per-month
]

globals [
  trans-year-threshold


  termination-node
  TimeToDiagnosis
  ;;scale-free distribution parameters
  degree-dist
  degree-dist-Bin
  degree-dist-all
  degree-dist-Bin-all
  exposures-per-month

  ;;track susceptibles not agents in an array of degree distribution
  susceptible-degree-dist

  num-exposed-sus ;Does not influence results, Maximum number of agents that could get infected at a timestep
  num-new-contacts ;Could be modified, Number of contacts needed to evolve an agent's network
  desired-neighbor-degree ;Could be modified, The number of contacts a neighbor of a newly infected agent should have
  SF-Distribution ;Scale-free degree-degree distribution
  avg-inf-degree ;Does not influence results, Average degree of newly infected agents
  newly-infected-degree-list ;For DD correlation data collection, could be deleted
  susceptible-contacts-degree-list ;For DD correlation data collection, could be deleted
  conditional-degree-dist-list ;For DD correlation data collection, could be deleted
  w  ;For NN
  xNN  ;For NN
  b1 ;For NN
  b2 ;For NN
  min_scale ;For NN
  max_scale ;For NN
  hidden_layer
  input
  max_degree
  avg_degree
  global_clustering_coefficient
  proportion_infected
  check_output

]

to testECNA
  file-close-all
  ; carefully [file-delete "TransmissionClusters.csv"] []
  ; write-cluster-statistics-header
  ;  carefully [file-delete "TransmissionClustersShort.csv"] []
  ;  write-cluster-header-Short
  set TimeToDiagnosis 10; upper bound for diagnosis in years
  setup-ECNA-globals
  setupECNA
  ;runECNA
  ; reportECNA
  ;  runClusterCodes
end
to runECNA
  let found false

  if count non-agent-susceptibles < 1000 [
    create-non-agent-susceptibles (termination-node * 2)[
    ]
  ]

  ;  ;; estimating ticks equivalent of number of infections in 10 years (say 2010 and 2015) assuming 50000 new infections per year
  ;  ;; trans-year will then use persons infected  in 10 years
  ;  if (found = false and count infecteds > (1 - (50000 * 10) / 1100000) * termination-node)[
  ;    set trans-year-threshold floor ((ticks + 1) / time-unit) / 2
  ;    set found true
  ;  ]
  goECNA
  repeat 10 [layout-ECNA]
  ;]

end

to setup-ECNA-globals

  set termination-node maxDegree * initial-infected * 10 ;500
  set exposures-per-month 100 / 12


  set degree-dist-all matrix:from-row-list [[0 0.214958105	0.154960174	0.213820213	0.175235337	0.126409434	0.072204407	0.030412744	0.011999586]]
  set degree-dist-Bin-all [0 1	2	4	8	16	32	64	128]; degree in each bin
  set degree-dist matrix:make-constant 1 length(degree-dist-Bin-all) 0
  normalizeDegreeDist 0


  set susceptible-degree-dist matrix:times-scalar degree-dist pop-size ;; determine number of persons in each bin
  let i 0
  while [i < length degree-dist-Bin][
    matrix:set susceptible-degree-dist 0 i round (matrix:get susceptible-degree-dist 0 i)
    set i i + 1
  ]

end

to clear
  clear-all
  file-close-all
end

to normalizeDegreeDist [sexInd]
  ;;NORMALIZE DEGREE DIST TO MAX DEGREE BINS (make it sum to 1)
  matrix:set-row degree-dist sexInd matrix:get-row degree-dist-all sexInd
  set degree-dist-Bin degree-dist-Bin-all
  let totalBins round(ln maxDegree / ln 2)

  let i  totalBins + 2
  repeat (length(degree-dist-Bin-all) - totalBins - 2)[
    matrix:set degree-dist sexInd i 0
    set degree-dist-Bin remove-item (length(degree-dist-Bin) - 1)  degree-dist-Bin
    set i i + 1
  ]

  set i 0
  let denominator 0
  repeat totalBins + 2[
    set denominator denominator + matrix:get degree-dist sexInd i
    set i i + 1
  ]
  let temp-row [] ;list
  set temp-row lput matrix:get-row degree-dist sexInd temp-row ;nested-list
  set temp-row matrix:from-row-list temp-row ;matrix
  set temp-row matrix:times-scalar temp-row (1 / denominator) ;matrix
  matrix:set-row degree-dist sexInd matrix:get-row temp-row 0

end
to get-mats
  ;let Scale-Free-Distribution csv:from-file "simulated_scale_free_dist_m5.csv" ;Theoretical scale-free distribution
  ;set SF-Distribution matrix:from-row-list Scale-Free-Distribution ;Distribution as a matrix

  let w1 csv:from-file "/Neural Network Weights/w_pref_attach_l.csv"
  set w matrix:from-row-list w1

  set xNN csv:from-file "/Neural Network Weights/x_pref_attach_l.csv"
  set xNN matrix:from-row-list xNN

  set b1 csv:from-file "/Neural Network Weights/b1_pref_attach_l.csv"
  set b1 matrix:from-row-list b1

  set b2 csv:from-file "/Neural Network Weights/b2_pref_attach_l.csv"
  set b2 matrix:from-row-list b2

  set min_scale csv:from-file "/Neural Network Weights/min_pref_attach_l.csv"
  set min_scale matrix:from-row-list min_scale

  set max_scale csv:from-file "/Neural Network Weights/max_pref_attach_l.csv"
  set max_scale matrix:from-row-list max_scale

  set hidden_layer item 0 matrix:dimensions xNN
  set input item 0 matrix:dimensions w
end

to-report return-random-degree
  let numAgents sum matrix:get-row susceptible-degree-dist 0
  let rand-degree random numAgents + 1
  let i 1
  let found false
  let cumulative 0
  while [i < (length degree-dist-Bin ) and found = false][
    set cumulative cumulative + matrix:get susceptible-degree-dist 0 i
    ; print cumulative
    if rand-degree <= cumulative [
      report item (i - 1) degree-dist-Bin + 1 + random (item (i) degree-dist-Bin - item (i - 1) degree-dist-Bin)
      set found true
      matrix:set susceptible-degree-dist 0 i ((matrix:get susceptible-degree-dist 0 i) - 1)
    ]
    set i i + 1
    ;  print i
  ]
end
to-report return-chosen-Bin [chosen-degree]
  let numAgents sum matrix:get-row susceptible-degree-dist 0
  let rand-degree random numAgents + 1

  let bin 0
  let i 0
  let found false
  while [i < length degree-dist-Bin and found = false][
    if chosen-degree <= item i degree-dist-Bin [
      ; item i degree-dist-Bin + 1 + random (item (i + 1) degree-dist-Bin - item i degree-dist-Bin )
      set bin i
      set found true
    ]
    set i i + 1
  ]

  set found false
  let chosen-bin 0
  ifelse ( matrix:get susceptible-degree-dist 0 bin + count-susceptible-agents bin) < 1[
    set i 1
    let chosenInd 0

    while [found = false and i < length degree-dist-Bin][
      let lower 0
      let upper 0
      carefully [set lower (matrix:get susceptible-degree-dist 0 bin - i) + count-susceptible-agents (bin - i)][set lower 0]
      carefully [set upper (matrix:get susceptible-degree-dist 0 bin + i) + count-susceptible-agents (bin + i)][set upper 0]
      ifelse random (lower + upper) + 1 <= lower [
        set chosenInd bin - i
      ]
      [set chosenInd bin + i
      ]

      if (lower != 0 or upper != 0)[
        set found true
        set chosen-bin chosenInd
      ]
      set i i + 1
    ]
  ]
  [
    set chosen-bin bin
    set found true
  ]

  report chosen-bin

end
to-report count-susceptible-agents [bin]
  let value 0
  carefully [set value count susceptibles with [desired-degree > item (bin - 1) degree-dist-Bin and desired-degree <= item bin degree-dist-Bin and (desired-degree - count my-links) > 0]]
  [set value 0]
  report value
end

to setupECNA
  reset-ticks
  carefully [file-delete "ECNAdata.csv"][]
  get-mats ;Stores all matrices for NN
  create-non-agent-susceptibles (termination-node * 10)[
  ]

  file-close

  create-infecteds initial-infected[
    set color red
    set size 0.5
    set shape "circle"
    set time-infected -1
    set aware? false
    set time-discovered random 100
    set desired-degree return-random-degree
    ;  set desired-degree random (item random  (item random length degree-list degree-list);file-read
    setxy random-xcor random-ycor
    set partners-per-month round (desired-degree / 12)
    if partners-per-month < 1 [set partners-per-month 1]
    ;print desired-degree
    ; print partners-per-month
    ;determine-non-eligible-modified
  ]

  ; kill-not-needed-links ;Kills all links a agent with no infected contacts has
  check-degree
  set max_degree last degree-dist-Bin

end

to check-degree ;Calculates the degree of each individual turtle
  ask infecteds[
    set degree count (my-links)
  ]
  ask susceptibles[
    set degree count (my-links)
  ]
end

to set-desired-degree
  ask infecteds [
    set desired-degree degree
  ]
  ask susceptibles [
    set desired-degree degree
  ]
end

to check-cluster ;Calculates the clustering-coefficient of each individual turtle
  ask infecteds[
    let clustering-coef nw:clustering-coefficient
    set clustering clustering-coef
  ]

  ask susceptibles[
    let clustering-coef nw:clustering-coefficient
    set clustering clustering-coef
  ]
end

to-report binomial [n p] ;Give a number of trials (n) and a probability (p) will compute an integer using a binomial distribution
  let bin_ct 0 ;Initializing counting variable
  repeat n [
    if random-float 1 < p
      [set bin_ct bin_ct + 1]
  ]
  report bin_ct
end
to writeECNA

  ask turtles with [breed = susceptibles][die]

  let n 0
  ask turtles with [ breed = susceptibles or breed = infecteds][
    set ID n
    set n n + 1
  ]

  file-open "ECNAdata.csv"
  ask infecteds[
    file-write who;ID
    file-write time-of-infection
    file-write time-of-diagnosis
    let source 0
    let source-inf-ID source-of-infection
    ask infecteds with [who = source-inf-ID]
    [set source ID]
    file-write source
    file-print ""
  ]
  file-close

end

to goECNA
  tick

  ;if (count infecteds >= (initial-infected + initial-susceptibles - 1) or num-exposed-sus = 0) [stop]
  ;if (ticks >= termination-ticks) [

  ask susceptibles [set contacts? false]
  ; determine-num-exposures ;Determines the number of exposures per time step for all agents
  infect-population-modified ;Determines if a susceptible agent will become infected
  calc-prop-infected
  determine-non-eligible-modified ;Determines which susceptible agents are not eligible to be linked with newly infected agents
  check-degree ;Calculates each turltes degree
               ;check-cluster ;Calculates each turtles clustering coefficient
  calc-avg-inf-degree ;Calculates average degree of newly infected agents
                      ; age-reset ;Updates turtles-own variables and global variables
                      ;reset ;Updates values
  check-awareness ;Checked the healthcare statues of all infected agents


  set num-exposed-sus 0
  ask infecteds[
    let sus-link link-neighbors with [breed = susceptibles]
    ask sus-link
    [set contacts? true]
  ]
  set num-exposed-sus count susceptibles with [contacts? = true]

end


to infect-population-modified

  ask infecteds with [count link-neighbors with [breed = susceptibles] > 0][
    let infected-node who
    let receiver n-of partners-per-month link-neighbors;Selects a random contact of aware-distributor
    let exposures-per-partner exposures-per-month / partners-per-month
    if (receiver != nobody ) [
      ask receiver[
        let prob-of-infection 1 - (1 - transmission-rate) ^ (exposures-per-partner) ;Calculates susceptible agents probability of infection
        let random-float1 random-float 1
        if (prob-of-infection > random-float1 and breed = susceptibles) [ ;Infects susceptible agents using the binomial distribution
          set breed  infecteds
          set time-of-infection ticks
          set source-of-infection infected-node
          let rand random-float 1

          ifelse rand < 0.25[;<0.7 years
            set time-of-diagnosis round(ticks + random (0.7 * time-unit));TimeToDiagnosis
          ]
          [ifelse rand < 0.5[;0.7 to 3 years
            set time-of-diagnosis round(ticks + 0.7 * time-unit + random (3 * time-unit - 0.7 * time-unit)); TimeToDiagnosis
            ]
            [
              ifelse rand < 0.75[; 3 to 7.8 years
                set time-of-diagnosis round(ticks + 3 * time-unit + random (7.8 * time-unit - 3 * time-unit));TimeToDiagnosis
              ]
              [;7.8 to TimeToDiagnosis
                set time-of-diagnosis round(ticks + 7.8 * time-unit + random (TimeToDiagnosis * time-unit - 7.8 * time-unit));
              ]

            ]
          ]


          set color red
          set size 0.5
          set shape "circle"
          set time-infected -1
          set aware? false
          set time-discovered random 100


          set partners-per-month round (desired-degree / 12)
          if partners-per-month < 1 [set partners-per-month 1]


        ]
      ]
    ]
  ]

end

to calc-prop-infected
  let num-inf count infecteds
  set proportion_infected (num-inf / pop-size) * 100
end


to determine-non-eligible-modified

  ask infecteds with [time-infected = -1 and desired-degree > 0] [ ;Susceptibles of the infected contacts of the newly infected agents are not eligible contacts
    set num-new-contacts 0

    ;pull from prob distribution (CDF) to determine degree
    let degree-of-newly-infected desired-degree
    let my_output calculate-dist-NN degree-of-newly-infected max_degree lambda avg_degree proportion_infected max_scale min_scale ;global_clustering_coefficient
    set check_output my_output
    let low-bound 0
    let high-bound matrix:get my_output 0 0

    set num-new-contacts desired-degree - count (my-links) ;The desired degree minus the current degree to determine number of links needed
                                                           ;print num-new-contacts
    let j 0
    while [j < num-new-contacts][
      let new-contacts nobody
      let desired-neighbor-degree1 desired-neighbor-degree-calc my_output low-bound high-bound
      let chosen-bin return-chosen-bin desired-neighbor-degree1

      if chosen-bin > 0[
        let eligible-agents  count-susceptible-agents (chosen-bin)
        let eligible-nonagents matrix:get susceptible-degree-dist 0 chosen-bin


        if eligible-agents + eligible-nonagents > 0[
          ifelse random (eligible-agents + eligible-nonagents) + 1 <= eligible-agents[
            set new-contacts one-of susceptibles with [desired-degree > item (chosen-bin - 1) degree-dist-Bin and desired-degree <= item chosen-bin degree-dist-Bin and (desired-degree - count my-links) > 0]
          ]
          [
            matrix:set susceptible-degree-dist 0 chosen-bin ((matrix:get susceptible-degree-dist 0 chosen-bin) - 1)
            set new-contacts one-of non-agent-susceptibles
            ask new-contacts[
              set breed susceptibles
              set desired-degree item (chosen-bin - 1) degree-dist-Bin + 1 + random (item (chosen-bin) degree-dist-Bin - item (chosen-bin - 1) degree-dist-Bin)
              set partners-per-month round (desired-degree / 12)
              if partners-per-month < 1 [set partners-per-month 1]

              set color green
              set size 0.5
              set time-in-simulation 0
            ]

          ]

        ]
      ]


      if (new-contacts != nobody) [
        create-ECNA-link-with new-contacts
      ]
      set j j + 1
    ]
  ]


end

to-report desired-neighbor-degree-calc [output_vector lb hb]
  let i 0
  let val random-float 1
  set desired-neighbor-degree 0
  while [ i < max_degree - 1] [ ;Draw from the theoretical distribution to determine the degree of the newly infected agent
    ifelse (val > lb) and (val < hb) [
      set desired-neighbor-degree i + 1
      set i i + 1
      set lb matrix:get output_vector (i - 1) (0)
      set hb matrix:get output_vector i (0)]

    [set i i + 1
      set lb matrix:get output_vector (i - 1) (0)
      set hb matrix:get output_vector i (0)]
  ]
  report desired-neighbor-degree
end

to-report calculate-dist-NN [node_degree maximum_degree network_lambda network_degree prop_inf maximum_scale minimum_scale] ; cluster_coef
  let output matrix:make-constant 1 max_degree matrix:get b2 0 0

  let j 0
  while [j < maximum_degree][       ;run for every neighbor degree (jth col)
    let h 0
    let mu (list node_degree (j + 1) network_lambda network_degree prop_inf) ;cluster_coef
                                                                             ;scale mu using (Input-min)/(max-min)
    let k 0
    while [k < length(mu)][
      set mu replace-item k mu ((item k mu - matrix:get minimum_scale 0 k) / (matrix:get maximum_scale 0 k - matrix:get minimum_scale 0 k))
      set k k + 1
    ]

    let v*h matrix:copy b1

    let v_h matrix:make-constant 1 hidden_layer 0

    while [h < hidden_layer][ ;run for every hidden network node (h=1:8)
                              ;summation of w(i,h)*mu(i) for every input (i=1:6) and add bias b1(h) [v*h]
                              ;mu(i) is the variable inpt corresponding to hidden layer node i
                              ;bias 1, b1(h) corresponds to every hidden layer node
      let i 0
      while [i < input][
        matrix:set v*h h 0 (matrix:get v*h h 0 + matrix:get w i h * item i mu)
        set i i + 1
      ]

      ;Apply activation function [v(h)]
      ;1/(1+e^-v*h)
      if matrix:get v*h h 0 <= 0[
        carefully [matrix:set v_h 0 h (1 / (1 + e ^ (- matrix:get v*h h 0)))] [matrix:set v_h 0 h 0]]
      if matrix:get v*h h 0 > 0[
        carefully [matrix:set v_h 0 h (1 / (1 + e ^ (- matrix:get v*h h 0)))] [matrix:set v_h 0 h 1]]
      ;summation of x(h)v(h) for every value of h (h=1:8) and add bias b2 [o]
      ;bias 2, b2, corresponds to the output
      ;store output in the jth col of vecotr

      matrix:set output 0 j  matrix:get output 0 j + matrix:get xNN h 0 * matrix:get v_h 0 h

      set h h + 1
    ]

    set j j + 1
  ]

  ;unscale output using output*(max-min)+min

  let k 0
  while [k < maximum_degree][
    matrix:set output 0 k (matrix:get output 0 k * (matrix:get maximum_scale 0 5 - matrix:get minimum_scale 0 5) + matrix:get minimum_scale 0 5)
    set k k + 1
  ]


  ; (normalize probabilities)
  let d 0
  let stopper 0
  let entry 0
  while [d < maximum_degree][
    set entry matrix:get output 0 d
    if entry < 0 [
      matrix:set output 0 d 0.00001]
    set d d + 1
  ]


  let row_sum 0
  let c 0
  while [c < maximum_degree][
    set row_sum row_sum + matrix:get output 0 c
    ;matrix:set output 0 c matrix:get output 0 c / sum item 0 matrix:to-row-list output
    set c c + 1
  ]
  set output output matrix:* (1 / row_sum)


  let m 0
  let mat_sum 0
  let output_cdf matrix:make-constant 1 max_degree 0

  while [m < maximum_degree][
    set mat_sum mat_sum + matrix:get output 0 m
    matrix:set output_cdf 0 m mat_sum
    set m m + 1
  ]


  set output_cdf matrix:transpose output_cdf
  report output_cdf
  ;report output
end

to calc-avg-inf-degree ;Calculates the average newly infected degree, does not influence model results
  set avg-inf-degree 0
  let total-links 0
  ask infecteds with [time-infected = -1][
    set total-links total-links + degree
  ]
  if (count infecteds with [time-infected = -1]) = 0 [stop]
  set avg-inf-degree (total-links / (count infecteds with [time-infected = -1]))
end

to-report avg-inf-degree1 ;Reports the average degree of newly infected agents, does not influence model results
  let avg-degree-val avg-inf-degree
  report avg-degree-val
end


to check-awareness ;Update healthcare status
  ask infecteds [
    if aware? = false and time-of-diagnosis > ticks[
      set aware? true
    ]
  ]
end

to kill-not-needed-links ;Kills links attached to agents with no infected contacts
  ask susceptibles[
    let my-contacts link-neighbors
    let my-contacts-links my-links
    if all? my-contacts [breed = susceptibles] [
      ask my-contacts-links [
        die
      ]
    ]
  ]

end

to-report global-clustering-coefficient ;Calculate global clustering coefficient
  let closed-triplets sum [ nw:clustering-coefficient * count my-links * (count my-links - 1) ] of turtles
  let triplets sum [ count my-links * (count my-links - 1) ] of turtles
  report closed-triplets / triplets
end

to-report infected-clustering-coefficient
  let sum-cc 0
  ask infecteds[
    set sum-cc sum-cc + clustering
  ]

  report sum-cc / (count infecteds)
end

to-report global-degree ;Calculate average global degree

  let total-links sum [degree] of infecteds
  let total-turtles count infecteds

  report (total-links) / total-turtles
end

to-report network-avg-degree
  let total-links sum [degree] of turtles
  let total-turtles count turtles
  report total-links / total-turtles
end



to-report median-degree
  let med-degree median [degree] of infecteds
  report med-degree
end

to-report variance-degree
  let var-degree variance [degree] of infecteds
  report var-degree
end


to-report median-degree-susc
  let med-degree median [count  link-neighbors with [breed = infecteds]] of susceptibles with [count link-neighbors with [breed = infecteds] > 0]
  report med-degree
end

to-report variance-degree-susc
  let var-degree variance [count  link-neighbors with [breed = infecteds]] of susceptibles with [count link-neighbors with [breed = infecteds] > 0]
  report var-degree
end

to layout-ECNA
  ; Some computational aspects of this layout function adopted from Netlogo Model Library codes:
  ; Wilensky, U. (1999). NetLogo. http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

  let agentset (turtle-set infecteds susceptibles)
  ;let agentset (turtle-set turtles)
  ask susceptibles [ifelse count my-links < 1 [set color black][set color green] set shape "circle" set size 0.5]
  let breed-links (link-set )
  ask infecteds [set breed-links (link-set breed-links my-links)]
  repeat 3 [
    ;; the more turtles we have to fit into the same amount of space,
    ;; the smaller the inputs to layout-spring we'll need to use
    let factor1 sqrt count agentset
    ;; numbers here are arbitrarily chosen for pleasing appearance
    layout-spring agentset breed-links (1 / factor1) (7 / factor1) (1 / factor1)
    display  ;; for smooth animation
  ]
  ;; don't bump the edges of the world
  let x-offset max [xcor] of agentset + min [xcor] of agentset
  let y-offset max [ycor] of agentset + min [ycor] of agentset
  ;; big jumps look funny, so only adjust a little each time
  set x-offset limit-magnitude1 x-offset 0.01
  set y-offset limit-magnitude1 y-offset 0.01
  ask agentset [ setxy (xcor - x-offset / 2) (ycor - y-offset / 2) ]
  ;if count Ptree-links > 0 [ ask Ttree-links [set thickness 0.5]]
  ;  if breed-links = Ptree-links [ ask Ptree-links [set thickness 0.2] ask Ttree-links [set thickness 1]]
end
to-report limit-magnitude1 [number limit]

  if number > limit [ report limit ]
  if number < (- limit) [ report (- limit) ]
  report number

end
@#$#@#$#@
GRAPHICS-WINDOW
375
17
1093
736
-1
-1
21.52
1
10
1
1
1
0
1
1
1
-16
16
-16
16
0
0
1
ticks
30.0

BUTTON
1
10
64
43
Clear
clear-all
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
70
10
138
43
setup
testECNA
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
6
429
95
474
NIL
count turtles
17
1
11

MONITOR
6
324
104
369
NIL
count infecteds
5
1
11

MONITOR
5
375
105
420
NIL
global-degree
0
1
11

MONITOR
108
326
232
371
SusceptibleAgents
count susceptibles with [count link-neighbors with [breed = infecteds] > 0]
0
1
11

PLOT
1025
23
1716
380
plot 1
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Num Infected" 1.0 0 -16777216 true "" "plot count infecteds"
"Num Susceptible Agents" 1.0 0 -7171555 true "" "plot num-exposed-sus"

PLOT
1029
392
1717
733
plot 2
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Avg Degree- of Infecteds" 1.0 0 -16777216 true "" "plot avg-inf-degree"

INPUTBOX
7
137
161
197
pop-size
1.0E7
1
0
Number

INPUTBOX
161
80
316
140
lambda
3.0
1
0
Number

INPUTBOX
161
139
312
199
initial-infected
5.0
1
0
Number

INPUTBOX
162
198
312
258
time-unit
12.0
1
0
Number

INPUTBOX
6
196
161
256
transmission-rate
0.01
1
0
Number

BUTTON
234
10
336
43
NIL
layout-ECNA
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
145
10
231
43
NIL
runECNA
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
5
78
160
138
maxDegree
8.0
1
0
Number

MONITOR
107
375
237
420
global-cluster-coefficient
global-clustering-coefficient
2
1
11

MONITOR
103
433
226
478
proportion_infected
proportion_infected
5
1
11

@#$#@#$#@
## WHAT IS IT?

Demonstration of agent-based evolving network modeling (ABENM) - a new simulation technique that combines features of agent-based modeling with compartmental modeling
It uses a new evolving contact network algorithm (ECNA) - A new algorithm for generation of scale-free networks that follow power-law degree distribution .

It is developed for simulation of diseases with low prevalence, where using current agent-based network modelingcan be computationaly infeasible, and where contact network structures are relevant to model that compartmnetal modeling is not suitable. 



## HOW TO USE IT
A basic understanidng of scale-free networks would be necessary to understand the concepts of this work. 

To run the model, on the interface, click on the following button, in order 
1. clear
2. setup
3. runECNA

## THINGS TO NOTICE

ABENM simulates only infected persons and their immediate contacts (infected and suceptible) as agents. All other susceptibles are simulated as compartmental modeling, with each compartment representing the degree (log-2 bin) of the nodes. 
As new persons becomes infected, the ECNA determines who the immediate contacts are- i.e., what is the their degree, and transitions correpsonding susceptible persons from the compartment to the network as neighbors of the newly infected node. 


## CREDITS AND REFERENCES

Developed by the disease modeling lab for https://blogs.umass.edu/chaitrag/projects/ 
Funding: NIH R01AI127236
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="ECNA" repetitions="100" runMetricsEveryStep="true">
    <setup>clear
setup</setup>
    <go>go</go>
    <metric>count infecteds</metric>
    <metric>global-degree</metric>
    <enumeratedValueSet variable="initial-prop-infected">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-nodes">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transmission-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportion-distributive-aware">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportion-distributive-unaware">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportion-recep">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reduction-factor">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="weight-of-triads">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="termination-ticks">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lambda">
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ABNM" repetitions="100" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>clear
setup-scale-free</setup>
    <go>run-abnm</go>
    <metric>count infecteds</metric>
    <metric>global-degree</metric>
    <enumeratedValueSet variable="initial-prop-infected">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-nodes">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transmission-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportion-distributive-aware">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportion-distributive-unaware">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportion-recep">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reduction-factor">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="weight-of-triads">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="termination-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lambda">
      <value value="1.79"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="cond_degree_dist" repetitions="1" runMetricsEveryStep="true">
    <setup>setup-scale-free</setup>
    <go>run-abnm</go>
    <metric>newly-infected-degree-list</metric>
    <metric>susceptible-contacts-degree-list</metric>
    <enumeratedValueSet variable="initial-infected">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transmission-rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maximum-degree">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportion-distributive-aware">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="weight-of-triads">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportion-distributive-unaware">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportion-recep">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reduction-factor">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-susceptibles">
      <value value="973"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="termination-ticks">
      <value value="251"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
