/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package com.rtg.assembler;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.rtg.mode.SequenceType;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SdfWriter;
import com.rtg.util.Constants;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.diagnostic.Diagnostic;

/**
 */
public final class CorrectReads {
  private CorrectReads() { }

  static void correct(File input, File output, int kmerSize, int thresholdFlag) throws IOException {
    try (final ReadPairSource source = ReadPairSource.makeSource(input, LongRange.NONE)) {
      final DeBruijnGraphBuilder builder = new DeBruijnGraphBuilder(Collections.singletonList(source), kmerSize, 0, 1);
      final DeBruijnGraph graph = builder.mDeBruijnGraph;
      final int threshold = thresholdFlag == -1 ? builder.calculateGoodThreshold() : thresholdFlag;
      Diagnostic.info("Using hash threshold: " + threshold);
      final List<SdfWriter> outputs = createOutputs(output, source);
      try {
        List<byte[]> fragments;
        final CorrectingMutator mutator = new CorrectingMutator();
        long readId = 0;
        source.reset();
        while ((fragments = source.nextFragments()) != null) {

          for (int i = 0; i < fragments.size(); ++i) {
            final byte[] fragment = fragments.get(i);

            final CorrectingMutator.SequenceBases best = correctRead(mutator, new CorrectingMutator.BaseRead(fragment), kmerSize, threshold, graph);
            final byte[] result;
            if (best == null) {
              result = fragment;
            } else {
              result = CorrectingMutator.sequenceBasesToBytes(best);
            }
//            System.err.println(DnaUtils.bytesToSequenceIncCG(fragment) + " -> " + DnaUtils.bytesToSequenceIncCG(result));
            final SdfWriter writer = outputs.get(i);
            writer.startSequence(String.valueOf(readId));
            writer.write(result, null, result.length);
            writer.endSequence();
          }

          ++readId;
        }
      } finally {
        for (SdfWriter writer : outputs) {
          writer.close();
        }

      }
    }
  }

  private static CorrectingMutator.SequenceBases correctRead(CorrectingMutator cm, CorrectingMutator.SequenceBases readBases, int kmerSize, int threshold, DeBruijnGraph graph) {
    final CorrectingMutator.SequenceBases output;
    CorrectingMutator.SequenceBases best = null;
    final int badHash;
    if ((badHash = firstBadHash(readBases, kmerSize, 0, threshold, graph)) != -1) {
      for (CorrectingMutator.SequenceBases mutation : cm.getMutations(readBases, badHash, readBases.length())) {
        if (firstBadHash(mutation, kmerSize, 0, threshold, graph) == -1) {
          if (best == null) {
            best = mutation;
            //System.err.println(reader.currentName() + " " + DnaUtils.bytesToSequenceIncCG(read) + " -> " + DnaUtils.bytesToSequenceIncCG(best));
          } else {
            best = null;
            break;
          }
        }
      }
      if (best != null) {
        output = best;
      } else {
        output = null;
      }
    } else {
      output = readBases;
    }
    return output;
  }

  private static int firstBadHash(CorrectingMutator.SequenceBases read, int kmerSize, int start, int threshold, DeBruijnGraph graph) {
    final KmerFactory fact = ByteKmer.factory();
    final byte[] bytes = CorrectingMutator.sequenceBasesToBytes(read);
    for (int i = start; i < read.length() - kmerSize; ++i) {
      final Kmer k = fact.make(bytes, i, i + kmerSize);
      if (!graph.contains(k)) {
        return i;
      } else if (graph.frequency(k) < threshold) {
        return i;
      }
    }
    return -1;
  }

  private static List<SdfWriter> createOutputs(File output, ReadPairSource source) {
    final List<String> outputNames = new ArrayList<>();
    if (source.numberFragments() == 2) {
      outputNames.add("left");
      outputNames.add("right");
    } else {
      for (int i = 0; i < source.numberFragments(); ++i) {
        outputNames.add(String.valueOf(i));
      }
    }
    final List<SdfWriter> outputs = new ArrayList<>(source.numberFragments());
    for (int i = 0; i < source.numberFragments(); ++i) {
      final SdfWriter writer = new SdfWriter(new File(output, outputNames.get(i)), Constants.MAX_FILE_SIZE, PrereadType.UNKNOWN, false, true, true, SequenceType.DNA);
      outputs.add(writer);
    }
    return outputs;
  }
}
