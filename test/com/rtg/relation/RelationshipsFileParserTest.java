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
package com.rtg.relation;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import com.rtg.relation.Relationship.RelationshipType;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class RelationshipsFileParserTest extends TestCase {

  private static final String RELATIONS_ODD = ""
    + "original-derived genomea genomeb randprop_mcprop=hoho=yoyo_foobar" + StringUtils.LS;

  public void testCases() throws IOException {
    final File dir = FileUtils.createTempDir("test", "relationshipfile");
    try {
      final File relationFile = new File(dir, "relationshipfile.txt");
      FileUtils.stringToFile(RELATIONS_ODD, relationFile);
      final GenomeRelationships gnf = RelationshipsFileParser.loadFile(relationFile);
      final String[] res = gnf.genomes();
      Arrays.sort(res);
      final Relationship[] expected = {new Relationship("genomea", "genomeb", RelationshipType.ORIGINAL_DERIVED)};
      expected[0].setProperty("randprop_mcprop", "hoho=yoyo_foobar");
      assertTrue(Arrays.toString(gnf.genomes()), Arrays.equals(new String[] {"genomea", "genomeb"}, res));
      assertTrue("Actual: " + Arrays.toString(gnf.relationships("genomea")), Arrays.equals(expected, gnf.relationships("genomea")));
      assertTrue("Actual: " + Arrays.toString(gnf.relationships("genomeb")), Arrays.equals(expected, gnf.relationships("genomeb")));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
  private static final String RELATIONS_BAD_3 = "" + "original-de genomea genomeb" + StringUtils.LS;
  private static final String RELATIONS_BAD_5 = "" + "original-derived genomeb" + StringUtils.LS;

  public void testBad() throws IOException {
    Diagnostic.setLogStream();
    final File dir = FileUtils.createTempDir("test", "relationshipfile");
    try {
      final File relationFile = new File(dir, "relationshipfile.txt");
      FileUtils.stringToFile(RELATIONS_BAD_3, relationFile);
      try {
        RelationshipsFileParser.loadFile(relationFile);
        fail("should throw exception");
      } catch (final NoTalkbackSlimException e) {
        assertEquals("unrecognized relationship type: 'original-de'", e.getMessage());
      }
      FileUtils.stringToFile(RELATIONS_BAD_5, relationFile);
      try {
        RelationshipsFileParser.loadFile(relationFile);
        fail("should throw exception");
      } catch (final NoTalkbackSlimException e) {
        assertEquals("unrecognized line in relationships: 'original-derived genomeb'", e.getMessage());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }


}
