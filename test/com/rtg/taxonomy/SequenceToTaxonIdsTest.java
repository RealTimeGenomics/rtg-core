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

package com.rtg.taxonomy;

import java.io.IOException;
import java.io.StringReader;
import java.util.Map;

import junit.framework.TestCase;

/**
 */
public class SequenceToTaxonIdsTest extends TestCase {

  /**
   * Test method for {@link com.rtg.taxonomy.SequenceToTaxonIds#sequenceToIds(java.io.File)}.
   * @throws IOException
   */
  public void testSequenceToIds() throws IOException {
    final String example = "1000\tmysequence1\n5454\tmysequence2 has extra stuff going on\n";
    Map<String, Integer> result = SequenceToTaxonIds.sequenceToIds(new StringReader(example));
    assertTrue(result.size() == 2);
    assertTrue(result.containsKey("mysequence1"));
    assertEquals(1000, (int) result.get("mysequence1"));
    assertEquals(5454, (int) result.get("mysequence2"));
  }

  public void testStatics() {
    assertEquals("taxonomy_lookup.tsv", TaxonomyUtils.TAXONOMY_TO_SEQUENCE_FILE);
  }
}
