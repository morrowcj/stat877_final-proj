# Response variables
## Insect Presence (0/1)
```
Harmandia
Phyllocolpa
Petiole.Gall
Leaf.Edge.Mine
Blotch.Mine
Lombardy.Mine
Weevil.Mine
Blackmine
Cottonwood.Leaf.Mine
Casebearer.Moth
Leafhoppers
Green.Aphids
Smokey.Aphids
Ants
Pale.Green.Notodontid
Aspen.Leaf.Beetle
Green.Sawfly
Cotton.Scale
```

## Insect density = `<insect>/Min.per.Tree`
See above

# Predictor Covariates
## Tree Traits
`Sex.Genet` (sex)
`PlantingDate` (age)
`Volume` (proxy for total growth)
`ALA.all` (leaf area, standardized)
`SLA.all` (leaf area, standardized)
`BBreakDegDay.all` (phenology, bud break time)
`EFNMean.all` (average number of extra-floral nectaries on the tree's leaves)

### leaf chemistry 
`CTsum` (condensed tannins %)
`PGsum` (phenolic glycosides %)
`Npct.all`(Nitrogen %)
`Cpct.all` (Carbon %)


## environmental factors
### ID variables
`Genet`
`id` (gwas sample names)
`SerialNo` (individual tree ID)
`Unique.ID` (observation ID, == <SerialNo>.<survey.event>)
`Block` (block within plot)
`Row` (row within plot)
`Position` (column within block)
`Latitude` (source lattitude of genet)
`Longitude` (source longitude of genet)


### Time
`Survey.Year` (year of survey)
`survey.event` (which survey event)
`Date` (date on which tree was surveyed)

### weather
`thunderstorm.event` (thunderstorm)
`avg.temp_F` (temperature)
`rain.event` (rain)
`fog.event` (fog)

# Offset
`Min.per.Tree` (time each tree was surveyed for)

