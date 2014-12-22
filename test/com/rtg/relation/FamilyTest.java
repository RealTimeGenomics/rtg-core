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

import static com.rtg.util.StringUtils.LS;

import java.io.BufferedReader;
import java.io.StringReader;
import java.util.Set;

import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class FamilyTest extends TestCase {

  /** Relations for testing. */
  public static final String RELATIONS = ""
      + "genome father sex=male disease=" + false + StringUtils.LS
      + "genome mother disease=" + true + StringUtils.LS
      + "genome childa disease=" + false + StringUtils.LS
      + "genome childb disease=" + true + StringUtils.LS
      + "parent-child father childa" + StringUtils.LS
      + "parent-child father childb" + StringUtils.LS
      + "parent-child mother childa" + StringUtils.LS
      + "parent-child mother childb" + StringUtils.LS
      ;


  public void testIndexes() {
    assertTrue(Family.FATHER_INDEX != Family.MOTHER_INDEX);
    assertTrue(Family.FATHER_INDEX != Family.FIRST_CHILD_INDEX);
    assertTrue(Family.MOTHER_INDEX != Family.FIRST_CHILD_INDEX);
  }

  public void testClassic() throws Exception {
    final Family f = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(RELATIONS))));
    assertTrue(f.getMother().equals("mother"));
    assertFalse(f.isDiseased(0));
    assertTrue(f.isDiseased(1));

    assertFalse(f.isDiseased("father"));
    assertTrue(f.isDiseased("mother"));
    assertEquals(2, f.getChildren().length);
    assertFalse(f.isDiseased("no-such-person"));
    assertTrue(f.isOneParentDiseased());

    assertFalse(f.isDiseased(2));
    assertTrue(f.isDiseased(3));
  }



  public void testMissingParent() throws Exception {
    final String input = ""
        + "parent-child mother childa" + LS
        + "parent-child mother childb" + LS
        ;
    try {
      Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(input))));
      fail();
    } catch (final IllegalArgumentException e) {
      // Expected
      assertEquals("There are fewer than two parents specified", e.getMessage());
    }
  }

  public void testChildHasOnlyOneParent() throws Exception {
    final String input = ""
        + "parent-child mother childa" + LS
        + "parent-child father childa" + LS
        + "parent-child father childb" + LS
        ;
    try {
      Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(input))));
      fail();
    } catch (final IllegalArgumentException e) {
      // Expected
      assertEquals("Child sample: 'childb' has 1 parents", e.getMessage());
    }
  }

  public void testScrewedUp() throws Exception {
    final String input = ""
        + "parent-child mother childa" + LS
        + "parent-child father childa" + LS
        + "parent-child father mother" + LS
        ;
    try {
      Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(input))));
      fail();
    } catch (final IllegalArgumentException e) {
      // Expected
    }
  }

  public void testScrewedUp2() throws Exception {
    final String input = ""
        + "parent-child mother childa" + LS
        + "parent-child mother father" + LS
        + "parent-child father childa" + LS
        ;
    try {
      Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(input))));
      fail();
    } catch (final IllegalArgumentException e) {
      // Expected
    }
  }

  public void testDoubleMother() throws Exception {
    final String input = ""
        + "parent-child mother childa" + LS
        + "parent-child mother childa" + LS
        + "parent-child father childb" + LS
        + "parent-child mother childb" + LS
        ;
    try {
      Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(input))));
      fail();
    } catch (final IllegalArgumentException e) {
      // Expected
      assertEquals("Child sample: 'childa' has 1 parents", e.getMessage());
    }
  }
  public void testTripleParents() throws Exception {
    final String input = ""
        + "parent-child mother childa" + LS
        + "parent-child mother childa" + LS
        + "parent-child father childb" + LS
        + "parent-child mother childb" + LS
        + "parent-child stepmother childb" + LS
        ;
    try {
      Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(input))));
      fail();
    } catch (final IllegalArgumentException e) {
      // Expected
      assertEquals("There are more than two parents specified", e.getMessage());
    }
  }

  public void testNormal() throws Exception {
    final String rel = "genome male sex=male" + StringUtils.LS
        + "genome female sex=female" + StringUtils.LS
        + "genome son sex=male" + StringUtils.LS
        + "parent-child male son" + StringUtils.LS
        + "parent-child female son" + StringUtils.LS
        + "genome daughter sex=female" + StringUtils.LS
        + "parent-child male daughter" + StringUtils.LS
        + "parent-child female daughter" + StringUtils.LS;
    Family fam = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(rel))));
    assertEquals("male", fam.getFather());
    assertEquals("female", fam.getMother());
  }

  public void testMultiFamily() {
    GenomeRelationships pedigree = new GenomeRelationships();
    pedigree.addGenome("father", GenomeRelationships.SEX_MALE).setProperty(GenomeRelationships.DISEASE_PROPERTY, "true");
    pedigree.addGenome("mother", GenomeRelationships.SEX_FEMALE);
    pedigree.addGenome("child", GenomeRelationships.SEX_MALE).setProperty(GenomeRelationships.DISEASE_PROPERTY, "true");
    pedigree.addParentChild("father", "child");
    pedigree.addParentChild("mother", "child");
    pedigree.addGenome("mother2", GenomeRelationships.SEX_FEMALE);
    pedigree.addGenome("child2", GenomeRelationships.SEX_MALE).setProperty(GenomeRelationships.DISEASE_PROPERTY, "true");
    pedigree.addGenome("child3", GenomeRelationships.SEX_MALE).setProperty(GenomeRelationships.DISEASE_PROPERTY, "true");
    pedigree.addParentChild("father", "child2");
    pedigree.addParentChild("mother2", "child2");
    pedigree.addParentChild("father", "child3");
    pedigree.addParentChild("mother2", "child3");
    pedigree.addGenome("unrelated1", GenomeRelationships.SEX_FEMALE);
    pedigree.addGenome("fatherfather", GenomeRelationships.SEX_MALE);
    pedigree.addGenome("fathermother", GenomeRelationships.SEX_FEMALE);
    pedigree.addParentChild("fatherfather", "father");
    pedigree.addParentChild("fathermother", "father");

    Set<Family> families = Family.getFamilies(pedigree, false, null);
    assertEquals(3, families.size());
    Family[] afamilies = families.toArray(new Family[families.size()]);
    assertEquals("father", afamilies[0].getFather());
    assertEquals("mother", afamilies[0].getMother());
    assertEquals("child", afamilies[0].getChildren()[0]);

    assertEquals("father", afamilies[1].getFather());
    assertEquals("mother2", afamilies[1].getMother());
    assertEquals("child2", afamilies[1].getChildren()[0]);
    assertEquals("child3", afamilies[1].getChildren()[1]);

    assertEquals("fatherfather", afamilies[2].getFather());
    assertEquals("fathermother", afamilies[2].getMother());
    assertEquals("father", afamilies[2].getChildren()[0]);
  }
}
