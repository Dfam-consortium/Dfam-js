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

   
  var a2m = seedAln.toA2M().alignments;

});

// Procesing a multiple seed-set Stockholm file simply requires
// a pre-processor to break up individual seed-sets on the "//" 
// separator.
