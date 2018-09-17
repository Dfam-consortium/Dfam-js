# Dfam-js Library 

## Introduction

A javascript libray for the Dfam/Dfam_consensus projects.

 
### DfamSeedAlignment
    
  A seed alignment is a set of related biological sequences 
  ( DNA, RNA, Amino Acids ) which have been aligned with respect
  to each other accounting for subtitions, deletions and insertions.
  A seed alignment is therefore a form of a sequence multiple
  alignment.  In Dfam we use seed alignments as the source of data
  for modeling a DNA sequence family using consensus sequences or
  as a profile Hidden Markov Model (HMM).  

  A Dfam seed alignment also contains metadata relavant to the
  family and specific to the domain of Transposable Elements.
  The object includes a parser for the standard seed alignment
  file format ( Stockholm ) as well as methods for converting
  the seed alignment to other formats ( ie. a summary format
  for visualization, and A2M format for storage in a database ).

#### Usage Example

Node:

```javascript
var DfamSeedAlignment = require('../src/DfamSeedAlignment.js');
var fs = require('fs');

// Create an instance of the class
var seedAln = new DfamSeedAlignment();

// Read in an example single seed-set Stockholm file
fs.readFile('example1.stk', 'utf8', function (err,data) {
  if (err)
    return console.log("Could not open file example1.stk for reading: " + err);

  // Parse the Stockholm data into the object
  seedAln.parseStockholm(data);

  // Convert the data to the AlignmentSummaryViewer format 
  var out = seedAln.toAlignmentSummary();
  console.log("Summary Data: " + JSON.stringify(out) );

  // Calculate a consensus
  out = seedAln.toConsensus();
  console.log("Consensus: " + out );

});
```

HTML:

```html
<script src="dist/Dfam-js.min.js"></script>

<script>

  var seedAln = new DfamSeedAlignment();

  var data = "...";  // <- Stockholm data

  // Parse the Stockholm data into the object
  seedAln.parseStockholm(data);

  // Convert the data to the AlignmentSummaryViewer format 
  var out = seedAln.toAlignmentSummary();
  console.log("Summary Data: " + JSON.stringify(out) );
  
  // Calculate a consensus
  out = seedAln.toConsensus();
  console.log("Consensus: " + out );

</script>
```
