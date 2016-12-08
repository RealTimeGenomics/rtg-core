/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.simulation.reads;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.rtg.reader.PrereadNamesInterface;
import com.rtg.reader.SequencesReader;
import com.rtg.util.AutoAddMap;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * We want to simulate meta-genomic communities
 * The input to this class is a map from taxonomy id to abundance or fraction of DNA
 * There may be more than one sequence associated with a single taxonomy id so we have to split the distribution across
 * those sequences
 */
public class TaxonomyDistribution {
  enum DistributionType {
    ABUNDANCE,
    DNA_FRACTION
  }
  private final double[] mDistribution;

  private static class TaxonMap extends AutoAddMap<Integer, TaxonSequences> {
    @Override
    public TaxonSequences make() {
      return new TaxonSequences();
    }
  }

  /**
   * Build a distribution
   * @param taxonomyDist input stream specifying the relative DNA fraction of taxonomy ids
   * @param taxonLookup mapping from sequence name to taxonomy id
   * @param reader reader providing length and name information
   * @param type specifies how to interpret the loaded distribution
   * @throws IOException if either the input stream parsing or sequence reader fails
   */
  TaxonomyDistribution(InputStream taxonomyDist, Map<String, Integer> taxonLookup, SequencesReader reader, DistributionType type) throws IOException {
    final Map<Integer, Double> taxonmyDistribution = parseTaxonDistribution(taxonomyDist);
    final AutoAddMap<Integer, TaxonSequences> taxon = new TaxonMap();
    final PrereadNamesInterface names = reader.names();
    final long numSeq = reader.numberSequences();
    final int[] lengths = reader.sequenceLengths(0, numSeq);
    for (int sequenceId = 0; sequenceId < numSeq; ++sequenceId) {
      final String name = names.name(sequenceId);
      final Integer taxonId = taxonLookup.get(name);
      final Double dist = taxonmyDistribution.get(taxonId);
      if (dist != null) {
        final TaxonSequences sequences = taxon.getOrAdd(taxonId);
        if (sequences != null) {
          sequences.add(new SequenceInfo(sequenceId, lengths[sequenceId]));
          sequences.mDist = dist;
        }
      }
    }
    final Map<Integer, Double> sequenceDist = new HashMap<>();
    for (TaxonSequences sequences : taxon.values()) {
      final Map<Integer, Double> dist;
      switch (type) {
       case ABUNDANCE:
        dist = sequences.simulationAbundanceDist();
         break;
       default:
         dist = sequences.simulationDnaFractionDist();
         break;
      }
      sequenceDist.putAll(dist);
    }
    Diagnostic.developerLog("Taxonomy distribution:" + sequenceDist.toString());
    mDistribution = new double[(int) reader.numberSequences()];
    double sum = 0;
    for (int i = 0; i < mDistribution.length; ++i) {
      final Double val = sequenceDist.get(i);
      mDistribution[i] =  val == null ? 0 : val;
      sum += mDistribution[i];
    }
    for (int i = 0; i < mDistribution.length; ++i) {
      mDistribution[i] = mDistribution[i] / sum;
    }


  }
  public double[] getDistribution() {
    return mDistribution;
  }

  /**
   *
   * @param is stream to read taxonomy distribution
   * @return map of taxonomy id to fraction
   * @throws IOException if the stream falls over
   */
  static Map<Integer, Double> parseTaxonDistribution(final InputStream is) throws IOException {
    HashMap<Integer, Double> map = new HashMap<>();
    double sum = 0;
    try (BufferedReader r = new BufferedReader(new InputStreamReader(is))) {
      String line;
      while ((line = r.readLine()) != null) {
        if (line.length() > 0 && line.charAt(0) != '#') {
          final String[] parts = line.trim().split("\\s+");
          if (parts.length != 2) {
            throw new IOException("Malformed line: " + line);
          }
          try {
            final double p = Double.parseDouble(parts[0]);
            if (p < 0 || p > 1) {
              throw new IOException("Malformed line: " + line);
            }
            sum += p;
            final Integer taxonId = Integer.parseInt(parts[1]);
            if (map.containsKey(taxonId)) {
              throw new IOException("Duplicated key: " + line);
            }
            map.put(taxonId, p);
          } catch (final NumberFormatException e) {
            throw new IOException("Malformed line: " + line);
          }
        }
      }
    }
    if (Math.abs(sum - 1) > 0.00001) {
      Diagnostic.warning("Input distribution sums to: " + String.format("%1.5g", sum));
      final Map<Integer, Double> oldmap = map;
      map = new HashMap<>();
      for (Map.Entry<Integer, Double> entry : oldmap.entrySet()) {
        map.put(entry.getKey(), entry.getValue() / sum);
      }
    }

    return map;
  }
  private static final class SequenceInfo {
    final int mSequenceId;
    final int mLength;

    private SequenceInfo(int id, int length) {
      mSequenceId = id;
      mLength = length;
    }
  }
  private static class TaxonSequences {
    double mDist;
    List<SequenceInfo> mSequences = new ArrayList<>();
    int totalLength() {
      int length = 0;
      for (SequenceInfo s : mSequences) {
        length += s.mLength;
      }
      return length;
    }
    void add(SequenceInfo s) {
      mSequences.add(s);
    }

    /**
     * @return <code>mDist</code> distributed amongst the underlying sequences in proportion to their lengths
     */
    Map<Integer, Double> simulationDnaFractionDist() {
      final Map<Integer, Double> res = new HashMap<>();
      final int total = totalLength();
      for (SequenceInfo s : mSequences) {
        res.put(s.mSequenceId, mDist * s.mLength / total);
      }
      return res;
    }
    /**
     * @return <code>mDist</code> distributed
     */
    Map<Integer, Double> simulationAbundanceDist() {
      final Map<Integer, Double> res = new HashMap<>();
      for (SequenceInfo s : mSequences) {
        res.put(s.mSequenceId, mDist * s.mLength);
      }
      return res;
    }
  }

}
