/*
 * Copyright (c) 2015. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.taxonomy;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reader.AbstractSdfWriter;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.MultiMap;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Utilities to help with loading taxonomy information from an SDF
 */
@TestClass("com.rtg.taxonomy.TaxStatsCliTest")
public final class TaxonomyUtils {

  /** Default name for taxonomy file. */
  public static final String TAXONOMY_FILE = "taxonomy.tsv";

  /** Default name of taxonomy to sequence name lookup file. */
  public static final String TAXONOMY_TO_SEQUENCE_FILE = "taxonomy_lookup.tsv";

  private TaxonomyUtils() { }

  /**
   * Returns true if the supplied reader contains taxonomy information
   * @param reader the input SDF
   * @return true if taxonomy information is present
   * @throws NoTalkbackSlimException if only one of the required taxonomy files is present
   */
  public static boolean hasTaxonomyInfo(SequencesReader reader) {
    final File taxonFile = new File(reader.path(), TAXONOMY_FILE);
    final File mappingFile = new File(reader.path(), TAXONOMY_TO_SEQUENCE_FILE);
    if (taxonFile.exists() && mappingFile.exists()) {
      return true;
    } else if (taxonFile.exists() || mappingFile.exists()) {
      throw new NoTalkbackSlimException("Reference SDF does not contain both taxonomy and sequences lookup");
    } else {
      return false;
    }
  }

  /**
   * Load a taxonomy from the supplied reader.
   * @param reader the input SDF
   * @return the taxonomy
   * @throws IOException if there are problems reading the input
   */
  public static Taxonomy loadTaxonomy(SequencesReader reader) throws IOException {
    final Taxonomy tax = new Taxonomy();
    try (FileInputStream fis = new FileInputStream(new File(reader.path(), TaxonomyUtils.TAXONOMY_FILE))) {
      tax.read(fis);
    }
    return tax;
  }

  /**
   * Load a map from (short) sequence names to taxon ids from the supplied reader
   * @param reader the input SDF
   * @return the map
   * @throws IOException if there are problems reading the input
   */
  public static Map<String, Integer> loadTaxonomyMapping(SequencesReader reader) throws IOException {
    return SequenceToTaxonIds.sequenceToIds(new File(reader.path(), TAXONOMY_TO_SEQUENCE_FILE));
  }

  /**
   * Load a map from taxon id to sequence ids from the supplied reader. Note that this implementation is
   * not efficient, in that it has to load reader names and invert the regular name to taxon id map.
   * @param reader the input SDF
   * @return the map
   * @throws IOException if there are problems reading the input
   */
  public static MultiMap<Integer, Long> loadTaxonomyIdMapping(SequencesReader reader) throws IOException {
    final AbstractSdfWriter.SequenceNameHandler handler = new AbstractSdfWriter.SequenceNameHandler();
    final Map<String, Long> names = ReaderUtils.getSequenceNameMap(reader);
    final Map<String, Integer> sequenceLookupMap = loadTaxonomyMapping(reader);
    // Invert the map and convert target to IDs
    final MultiMap<Integer, Long> result = new MultiMap<>();
    for (Map.Entry<String, Integer> entry : sequenceLookupMap.entrySet()) {
      final Long id = names.get(handler.handleSequenceName(entry.getKey()).label());
      if (id != null) {
        result.put(entry.getValue(), id);
      }
    }
    return result;
  }
}
