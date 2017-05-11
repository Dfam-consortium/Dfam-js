// 
//  DfamSeedAlignment 
//    
//    A seed alignment is a set of related biological sequences 
//    ( DNA, RNA, Amino Acids ) which have been aligned with respect
//    to each other accounting for subtitions, deletions and insertions.
//    A seed alignment is therefore a form of a sequence multiple
//    alignment.  In Dfam we use seed alignments as the source of data
//    for modeling a DNA sequence family using consensus sequences or
//    as a profile Hidden Markov Model (HMM).  
//
//    A Dfam seed alignment also contains metadata relavant to the
//    family and specific to the domain of Transposable Elements.
//    The object includes a parser for the standard seed alignment
//    file format ( Stockholm ) as well as methods for converting
//    the seed alignment to other formats ( ie. a summary format
//    for visualization, and A2M format for storage in a database ).
//
// Robert Hubley 2017
//    
(function(global) {
  'use strict';

  function DfamSeedAlignment() {
    this.headers = {};
    this.alignments = [];
  }

  DfamSeedAlignment.prototype.parseStockholm = function (textData) {
        if (textData == null)
            return;
        var lines = textData.split(/\r\n|\n/);
        var stockholmObj = this;
        for (var i = 0; i < lines.length; i++) {
            var matches;
            if ((matches = lines[i].match(/^\s*$/)) != null) {
                // Blank lines
            } else if ((matches = lines[i].match(/^#\s+STOCKHOLM\s+([\.\d]+)/)) != null) {
                stockholmObj.version = matches[1];
            } else if ((matches = lines[i].match(/^#=GF\s+(\S+)\s+(\S+)/)) != null) {
                // Need more sophistication here for various types of tags. i.e
                // CC is an array type tag
                stockholmObj.headers[matches[1]] = matches[2];
            } else if ((matches = lines[i].match(/^#=GC\s+RF\s+(\S+)/)) != null) {
                stockholmObj.headers['RF'] = matches[1];
            } else if ((matches = lines[i].match(/^\s*(\S+)\s+([ACGTUMRWSYKVHDBNacgtumrwsykvhdbn\.]+)\s*$/)) != null) {
                stockholmObj.alignments.push({
                    seqid: matches[1],
                    alignment: matches[2]
                });
            } else if ((matches = lines[i].match(/^\/\/\s*$/)) != null) {
                // End of seed alignment
                break;
            } else {
                console.log('parseStockholmData: Unknown line ' + i + ': ' + lines[i]);
            }
        }
        if (stockholmObj.headers.RF == null)
            throw 'Stockholm Alignment: Missing RF line!  Not a consensus based multiple alignment';
        // console.log("OBJ: " + JSON.stringify(stockholmObj, null, '\t') );
        return stockholmObj;
    };

    // species:assembly:sequence:start-end
    // assembly:sequence:start-end
    // sequence:start-end
    // legacy:
    //    sequence_start_end_[F/R]  Not exactly what we want.  I.e the start/end
    //    do not refer to the aligned sequence range rather it's the contig that it
    //    was found in.  So if we want to see what "should" have aligned we should
    //    pass along more metadata.
    //
    //   Where orientation is encoded by the order of the start/end values.  I.e
    //   start < end  ( orientation = forward strand )
    //   end < start ( orientation = reverse strand )
    //   start = end ( orientation can only be determined by referencing the source seq )
    DfamSeedAlignment.prototype.toAlignmentSummary = function () {
        var stockholmObj = this;
        var consensus = this.generateConsensusFromAlignment2();
        var aligns = stockholmObj.alignments;
        var windowSize = 10;

        var summaryData = {};
        summaryData.qualityBlockLen = windowSize;
        summaryData.length = (consensus.match(/[^\.]/g) || []).length;
        summaryData.alignments = [];

        // For each alignment
        for (var i = 0; i < aligns.length; i++) {
            // Decode ID
            // species:assembly:sequence:start-end
            var idNomenclatureRE = /^(?:([^:\s]+):)?(?:([^:\s]+):)?(\S+):(\d+)-(\d+)$/;
            // sequence_start_end_[F/R]
            var legacyNomenclatureRE = /^(\S+)_(\d+)_(\d+)(?:_(R))?$/;
            var matches;
            var seqID;
            var seqStart = 0;
            var seqEnd = 0;
            var orient = 'F';
            // Don't include the reference sequence in the alignment.
            if (/^ref:/.test(aligns[i].seqid))
                continue;
            if ((matches = aligns[i].seqid.match(idNomenclatureRE)) != null) {
                seqID = matches[3];
                if (matches[4] <= matches[5]) {
                    seqStart = matches[4];
                    seqEnd = matches[5];
                    orient = 'F';
                } else {
                    seqStart = matches[5];
                    seqEnd = matches[4];
                    orient = 'R';
                }
            } else if ((matches = aligns[i].seqid.match(legacyNomenclatureRE)) != null) {
                // Nomenclature used by Linup
                seqID = matches[1];
                seqStart = matches[2];
                seqEnd = matches[3];
                orient = 'F';
                if (matches[4] == 'R')
                    orient = 'R';
            } else {
                // throw "generateAlignmentHeatMapData: Stockholm file is using a strange nomenclature. id = " + aligns[i].seqid;
                // Tolerance
                seqID = aligns[i].seqid;
                seqStart = 0;
                seqEnd = 0;
                orient = '';
            }
            var ins = 0;
            var del = 0;
            var aligned = 0;
            var mut = 0;
            var scores = [];
            var conStart = 0;
            // var conEnd = 0;
            var consPos = 0;
            var startMatch = aligns[i].alignment.match(/^(\.*)[^\.]/);
            var endMatch = aligns[i].alignment.match(/.*[^\.](\.*)$/);
            var alignStart = 0;
            if (startMatch.length > 0)
                alignStart = startMatch[1].length;
            var alignEnd;
            if (endMatch.length > 0)
                alignEnd = aligns[i].alignment.length - endMatch[1].length - 1;
            // For each column where consensus contains a base
            for (var j = 0; j <= alignEnd; j++) {
                var cBase = consensus.substring(j, j + 1);
                if (cBase != '.') {
                    consPos++;
                }

                if (j < alignStart)
                    continue;

                // Now we can set the consensus start position
                if (conStart == 0)
                    conStart = consPos;

                var aBase = aligns[i].alignment.substring(j, j + 1);
                if (cBase == '.') {
                    if (aBase != '.')
                        ins = ins + 1;
                    continue;
                }

                aligned = aligned + 1;
                if (aBase == '.') {
                    del = del + 1;
                } else {
                    if (conStart == 0)
                        conStart = aligned;
                    if (cBase != aBase)
                        mut = mut + 1;
                }
                if (aligned % windowSize == 0) {
                    var score = windowSize - (mut + del);
                    if (ins)
                        score = score - 1;
                    if (score < 1)
                        score = 1;
                    scores.push(score);
                    ins = 0;
                    del = 0;
                    mut = 0;
                }
            }
            if (mut || del || ins) {
                var score = windowSize - (mut + del);
                if (ins)
                    score = score - 1;
                if (score < 1)
                    score = 1;
                scores.push(score);
            }
            var record = [];
            record.push(seqID);
            record.push(conStart);
            record.push(aligned);
            record.push(scores);
            record.push(orient);
            record.push('0.0'); // TODO: divergence...
            record.push(seqStart); // TODO: Document the source of this
            record.push(seqEnd);
            summaryData.alignments.push(record);
        }
        summaryData.num_alignments = summaryData.alignments.length;
        return (summaryData);
    };

    // Majority rule simple consensus ( RF-only alignments, ie. consensus based )
    DfamSeedAlignment.prototype.toConsensus = function() {
        var stockholmObj = this;
        var rfData = stockholmObj.headers.RF;
        var aligns = stockholmObj.alignments;
        var colData = [];
        for (var i = 0; i < aligns.length; i++) {
            for (var j = 0; j < aligns[i].alignment.length; j++) {
                var aBase = aligns[i].alignment.substring(j, j + 1);
                if (colData[j] == null)
                    colData[j] = {};
                if (aBase != '.') {
                    if (colData[j][aBase] == null)
                        colData[j][aBase] = 0;
                    colData[j][aBase]++;
                }
            }
        }
        var consensus = '';
        for (var j = 0; j < rfData.length; j++) {
            if (rfData.substring(j, j + 1) == 'x') {
                var highCount = 0;
                var highBase = '';
                Object.keys(colData[j]).forEach(function(aBase) {
                    if (colData[j][aBase] > highCount) {
                        highCount = colData[j][aBase];
                        highBase = aBase;
                    }
                });
                consensus = consensus + highBase;
            }
        }
        return consensus;
    };

    //   Refine the consensus sequence given the multiple alignment data.
    //   Correct for missed CpG calls.  Note: The correction currently
    //   assumes a AT bias in the genome ( good for mammals ) in the
    //   hardcoded lineup matrix.
    DfamSeedAlignment.prototype.generateConsensusFromAlignment2 = function(stockholmObj) {
        var stockholmObj = this;
        //
        // This is what boosts the value of the CG
        // consensus comparison if the instance dinucs
        // are CA or TG.
        //
        var CGParam = 12; // TG or CA match is 19, TG <-> CA mismatch -12
        // Previously set at 14. Seems to overestimate in very
        //   old elements.
        var TAParam = -5; // TG or CA to TA is -4, so slightly worse than that
        //   CG -> TA mismatch would have been -8
        var CGTransParam = 2; // Adjust scores of transitions/transversion pairs
        //  that could have arisen from CpG site.

        // For mammals where there is a strong A/T bias
        //  A   R   G   C   Y   T   K   M   S   W   N   X   Z
        var matrix = [
            [9, 0, -8, -15, -16, -17, -13, -3, -11, -4, -2, -7, -3, -6],
            [2, 1, 1, -15, -15, -16, -7, -6, -6, -7, -2, -7, -3, -6],
            [-4, 3, 10, -14, -14, -15, -2, -9, -2, -9, -2, -7, -3, -6],
            [-15, -14, -14, 10, 3, -4, -9, -2, -2, -9, -2, -7, -3, -6],
            [-16, -15, -15, 1, 1, 2, -6, -7, -6, -7, -2, -7, -3, -6],
            [-17, -16, -15, -8, 0, 9, -3, -13, -11, -4, -2, -7, -3, -6],
            [-11, -6, -2, -11, -7, -3, -2, -11, -6, -7, -2, -7, -3, -6],
            [-3, -7, -11, -2, -6, -11, -11, -2, -6, -7, -2, -7, -3, -6],
            [-9, -5, -2, -2, -5, -9, -5, -5, -2, -9, -2, -7, -3, -6],
            [-4, -8, -11, -11, -8, -4, -8, -8, -11, -4, -2, -7, -3, -6],
            [-2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -1, -7, -3, -6],
            [-7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -3, -6],
            [-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -6],
            [-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, 3]
        ];
        var alphabet_r = ['A', 'R', 'G', 'C', 'Y', 'T', 'K', 'M', 'S', 'W', 'N', 'X', 'Z', '.'];
        var alphabet_h = {};
        for (var i = 0; i < alphabet_r.length; i++) {
            alphabet_h[alphabet_r[i]] = i;
        }

        // var rfData = stockholmObj.headers.RF;
        var aligns = stockholmObj.alignments;
        var colData = [];
        for (var i = 0; i < aligns.length; i++) {
            var startMatch = aligns[i].alignment.match(/^(\.*)[^\.]/);
            var endMatch = aligns[i].alignment.match(/.*[^\.](\.*)$/);
            var alignStart = 0;
            if (startMatch.length > 0)
                alignStart = startMatch[1].length;
            var alignEnd;
            if (endMatch.length > 0)
                alignEnd = aligns[i].alignment.length - endMatch[1].length - 1;
            for (var j = alignStart; j <= alignEnd; j++) {
                var aBase = aligns[i].alignment.substring(j, j + 1);
                if (colData[j] == null)
                    colData[j] = {};
                if (colData[j][aBase] == null)
                    colData[j][aBase] = 0;
                colData[j][aBase]++;
            }
        }

        // Calculate as-is score
        var posHighScores = [];
        var consensus = '';
        // Foreach column
        for (var i = 0; i < aligns[0].alignment.length; i++) {
            var maxScore = -1000000000;
            var highBase;
            if (colData[i] != null) {
                // Foreach possible consensus base
                for (var j = 0; j < alphabet_r.length; j++) {
                    var cons = alphabet_r[j];
                    var score = 0;
                    // Foreach observed base
                    Object.keys(colData[i]).forEach(function(obs) {
                        // score += count_obs * matrix( cons, obs )
                        score = score + (colData[i][obs] * matrix[j][alphabet_h[obs]]);
                    });
                    // TODO: optimisation...make sure stockholm parser upperCases
                    // all alignment data
                    if (score == maxScore && /[ACGTacgt]/.test(cons)) {
                        highBase = cons;
                    } else if (score > maxScore) {
                        highBase = cons;
                        maxScore = score;
                    }
                }
                consensus = consensus + highBase;
                posHighScores.push(maxScore);
            }
        }

        // for each pair of dinucleotides in the consensus
        var diRE = /([^\.])(\.*)/g;
        var matches;
        while ((matches = diRE.exec(consensus)) != null) {
            var rgtIdx = diRE.lastIndex;
            var lftIdx = rgtIdx - matches[2].length - 1;
            var consLeft = matches[1];
            var consRight = consensus.substring(rgtIdx, rgtIdx + 1);
            // var consDiNucl = consLeft + consRight;
            var CGScore = 0;
            var dnScore = 0;
            if (consRight != '') {
                // console.log("consDi = " + consDiNucl);
                for (var i = 0; i < aligns.length; i++) {
                    var hitLeft = aligns[i].alignment.substring(lftIdx, lftIdx + 1);
                    var hitRight = aligns[i].alignment.substring(rgtIdx, rgtIdx + 1);
                    if (hitLeft == '.' || hitRight == '.')
                        continue;
                    var hitDiNucl = hitLeft + hitRight;
                    // console.log("   hitDiNucl = " + hitDiNucl);
                    // Recalculate the score of the consensus ( excluding gap characters )
                    dnScore += matrix[alphabet_h[consLeft]][alphabet_h[hitLeft]];
                    dnScore += matrix[alphabet_h[consRight]][alphabet_h[hitRight]];
                    if (hitDiNucl == 'CA' || hitDiNucl == 'TG') {
                        CGScore = CGScore + CGParam;
                    } else if (hitDiNucl == 'TA') {
                        CGScore = CGScore + TAParam;
                    } else if (hitDiNucl == 'TC' || hitDiNucl == 'TT') {
                        // in other words; C->T transition scores +2
                        // transversion scored normally
                        CGScore = CGScore + CGTransParam + (matrix[alphabet_h['G']][alphabet_h[hitRight]]);
                    } else if (hitDiNucl == 'AA' || hitDiNucl == 'GA') {
                        CGScore = CGScore + CGTransParam + (matrix[alphabet_h['C']][alphabet_h[hitLeft]]);
                        // same as above
                    } else {
                        CGScore = CGScore + matrix[alphabet_h['C']][alphabet_h[hitLeft]];
                        CGScore = CGScore + matrix[alphabet_h['G']][alphabet_h[hitRight]];
                    }
                }
                if (CGScore > dnScore) {
                    consensus = consensus.substring(0, lftIdx) + 'C' + consensus.substring(lftIdx + 1, consensus.length);
                    consensus = consensus.substring(0, rgtIdx) + 'G' + consensus.substring(rgtIdx + 1, consensus.length);
                }
            }
        }
        return consensus;
    };

    //
    // A2M Format
    //
    //  From https://compbio.soe.ucsc.edu/a2m-desc.html:
    //
    //  "Uppercase characters and "-" represent alignment columns, and there must
    //   be exactly the same number of alignment columns per sequence.  Lowercase
    //   characters ( and spaces or "." ) represent insertion positions between
    //   alignment columns or at the ends of the sequence.  The spaces or periods
    //   in the multiple alignments are only for human readability, and
    //   may be omitted."
    //
    //   A3M has been used informally to designate files with the "." omitted.
    //
    // For example the following A2M:
    //
    // >seq1
    // AAACCTagtgCGGGATC
    // >seq2
    // ----CTCG-TATC
    // >seq3
    // --AC--CGGGATC
    // >seq4
    // AAAaCCTCG-----
    //
    // Represents the full stockholm multiple alignment:
    //  cons:   AAA.CCT....CGGGATC
    //  --------------------------
    //  seq1:   AAA.CCTAGTGCGGGATC
    //  seq2:   .....CT....CG.TATC
    //  seq3:   ..A.C......CGGGATC
    //  seq4:   AAAACCT....CG.....
    //
    DfamSeedAlignment.prototype.toA2M = function() {
        var stockholmObj = this;
        if (stockholmObj == null)
            return;

        var a2mObj = {
           alignments: []
        };

        // The RF line may already be provided in the stockholm
        // file.
        var RF;
        if (stockholmObj.headers != null &&
            stockholmObj.headers.RF != null)
            RF = stockholmObj.headers.RF;

        // If the RF line wasn't present we must infer the
        // insertions by generating a consensus for the
        // multiple alignment.
        if (RF == null) {
            var cons = this.generateConsensusFromAlignment2();
            RF = cons.replace(/[^\.]/g, 'X');
        }

        var aligns = stockholmObj.alignments;
        var matchColCnt = -1;
        for (var i = 0; i < aligns.length; i++) {
            // Create a copy of the source data
            var a2m_seq = aligns[i].alignment;

            // Uppercase all seqs and convert all "-"s to "."s to
            // have a consistent starting point.
            a2m_seq = a2m_seq.toUpperCase();
            a2m_seq = a2m_seq.replace(/-/g, '.');

            var seqLen = (a2m_seq.match(/[^\.]/g) || []).length;

            // Now go over each column in the sequence.
            var rfPos = 0;
            var modelStart = -1;
            var modelEnd = -1;
            for (var j = 0; j < a2m_seq.length; j++) {
                // Get sequence character at this position
                var base = a2m_seq.substring(j, j + 1);
                // Get RF data for this column.
                var rf = RF.substring(j, j + 1);
                if (rf == '.') {
                    if (base != '.') {
                        a2m_seq = a2m_seq.substr(0, j) +
                            base.toLowerCase() +
                            a2m_seq.substr(j + 1);
                    }
                } else {
                    rfPos++;
                    if (base == '.') {
                        a2m_seq = a2m_seq.substr(0, j) +
                            '-' +
                            a2m_seq.substr(j + 1);
                    }else {
                      if (modelStart < 0 )
                        modelStart = rfPos;
                      modelEnd = rfPos;
                    }
                }
            }
            // Remove remaining "."s
            a2m_seq = a2m_seq.replace(/\./g, '');
            var tmpCnt = (a2m_seq.match(/[^A-Z\-]/g) || []).length;
            if (matchColCnt >= 0 && matchColCnt != tmpCnt)
                throw 'stockholmToA2M: Error converting to A2M format.  The ' +
                    'number of match columns is not consistent ( was ' + matchColCnt +
                    ' and now ' + tmpCnt + ' ).  Here is the offending line:' +
                    aligns[i].alignment;
            a2mObj.alignments.push({
                seq_id: aligns[i].seqid,
                seq_start: aligns[i].start || 1,
                seq_end: aligns[i].end || seqLen,
                a2m_seq: a2m_seq,
                strand: aligns[i].orient || '+',
                model_start: modelStart,
                model_end: modelEnd
            });
        }
        return a2mObj;
    };


  if (typeof module === 'object' && module && typeof module.exports === 'object') {
    // Expose functions/objects for loaders that implement the Node module pattern.
    module.exports = DfamSeedAlignment;
  } else {
    // Otherwise expose ourselves directly to the global object.
    global.DfamSeedAlignment = DfamSeedAlignment;
    // Register as a named AMD module.
    if (typeof define === 'function' && define.amd) {
      define('dfamseedalignment', [], function() {
        return DfamSeedAlignment;
      });
    }
  }
}(this.window || (typeof global != 'undefined' && global) || this));
