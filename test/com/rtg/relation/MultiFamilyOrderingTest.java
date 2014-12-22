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

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import junit.framework.TestCase;

/**
 */
public class MultiFamilyOrderingTest extends TestCase {

  public void testSuperDuperUltraSimple() {
    final Set<Family> fams = new HashSet<>();
    final Family fam = new Family("father", "mother", "child1", "child2", "child3");
    fams.add(fam);
    final List<Family> famsO = MultiFamilyOrdering.orderFamiliesAndSetMates(fams);
    assertEquals(1, famsO.size());
    assertTrue(fams.contains(famsO.get(0)));
    assertEquals(1, fams.size());
    assertTrue(MultiFamilyOrdering.isMonogamous(famsO));
    assertEquals(1, fam.getFatherDistinctMates());
    assertEquals(1, fam.getMotherDistinctMates());
    //System.err.println(famsO);
  }

  public void testMultiGen() {
    final Set<Family> fams = new LinkedHashSet<>();
    final Family gen1 = new Family("father", "mother", "child1", "child2", "child3");
    final Family gen2 = new Family("child1", "dInLaw1", "gchild1", "gchild2");
    fams.add(gen2);
    fams.add(gen1);
    final List<Family> famsO = MultiFamilyOrdering.orderFamiliesAndSetMates(fams);
    assertEquals(2, famsO.size());
    assertEquals(gen1, famsO.get(0));
    assertEquals(gen2, famsO.get(1));
    fams.clear();
    fams.add(gen1);
    fams.add(gen2);
    assertEquals(2, famsO.size());
    assertEquals(gen1, famsO.get(0));
    assertEquals(gen2, famsO.get(1));
    assertTrue(MultiFamilyOrdering.isMonogamous(famsO));
    assertEquals(1, gen1.getFatherDistinctMates());
    assertEquals(1, gen1.getMotherDistinctMates());
    assertEquals(1, gen2.getFatherDistinctMates());
    assertEquals(1, gen2.getMotherDistinctMates());
  }

  public void testMultiGenExtra() {
    final Set<Family> fams = new LinkedHashSet<>();
    final Family gen1 = new Family("father", "mother", "child1", "child2", "child3");
    final Family gen2 = new Family("child1", "dInLaw1", "gchild1", "gchild2");
    final Family gen2b = new Family("sInLaw1", "child3", "gchild4", "gchild5");
    fams.add(gen2);
    fams.add(gen2b);
    fams.add(gen1);
    final List<Family> famsO = MultiFamilyOrdering.orderFamiliesAndSetMates(fams);
    assertEquals(3, famsO.size());
    assertEquals(gen1, famsO.get(0));
    assertTrue(famsO.contains(gen2));
    assertTrue(famsO.contains(gen2b));
    assertTrue(MultiFamilyOrdering.isMonogamous(famsO));
    assertEquals(1, gen1.getFatherDistinctMates());
    assertEquals(1, gen1.getMotherDistinctMates());
    assertEquals(1, gen2.getFatherDistinctMates());
    assertEquals(1, gen2.getMotherDistinctMates());
    assertEquals(1, gen2b.getFatherDistinctMates());
    assertEquals(1, gen2b.getMotherDistinctMates());
  }

  public void testMultiGenDifferentFamilies() {
    final Set<Family> fams = new LinkedHashSet<>();
    final Family gen1 = new Family("father", "mother", "child1", "child2", "child3");
    final Family gen2 = new Family("child1", "dInLaw1", "gchild1", "gchild2");
    final Family gen1b = new Family("barney", "betty", "bambam");
    final Family gen2b = new Family("bambam", "pebbles", "gchild4");
    fams.add(gen1);
    fams.add(gen2);
    fams.add(gen1b);
    fams.add(gen2b);
    final List<Family> famsO = MultiFamilyOrdering.orderFamiliesAndSetMates(fams);
    assertEquals(4, famsO.size());
    assertTrue(famsO.containsAll(Arrays.asList(gen1, gen1b, gen2, gen2b)));
    assertTrue(famsO.indexOf(gen1) < famsO.indexOf(gen2));
    assertTrue(famsO.indexOf(gen1b) < famsO.indexOf(gen2b));
    assertTrue(MultiFamilyOrdering.isMonogamous(famsO));
    assertEquals(1, gen1.getFatherDistinctMates());
    assertEquals(1, gen1.getMotherDistinctMates());
    assertEquals(1, gen2.getFatherDistinctMates());
    assertEquals(1, gen2.getMotherDistinctMates());
    assertEquals(1, gen1b.getFatherDistinctMates());
    assertEquals(1, gen1b.getMotherDistinctMates());
    assertEquals(1, gen2b.getFatherDistinctMates());
    assertEquals(1, gen2b.getMotherDistinctMates());
  }

  public void testGG() {
    final Set<Family> fams = new LinkedHashSet<>();
    final Family gen1 = new Family("father", "mother", "child1", "child2", "child3");
    final Family gen2 = new Family("child1", "dInLaw1", "gchild1", "gchild2");
    final Family gen1b = new Family("barney", "betty", "bambam");
    final Family gen2b = new Family("bambam", "pebbles", "gchild4");
    final Family gen3 = new Family("gchild4", "gchild2", "ggchild0");
    fams.add(gen1);
    fams.add(gen2);
    fams.add(gen3);
    fams.add(gen1b);
    fams.add(gen2b);
    final List<Family> famsO = MultiFamilyOrdering.orderFamiliesAndSetMates(fams);
    assertEquals(5, famsO.size());
    assertTrue(famsO.containsAll(Arrays.asList(gen1, gen1b, gen2, gen2b, gen3)));
    assertTrue(famsO.indexOf(gen1) < famsO.indexOf(gen2));
    assertTrue(famsO.indexOf(gen1b) < famsO.indexOf(gen2b));
    assertTrue(famsO.indexOf(gen2b) < famsO.indexOf(gen3));
    assertTrue(famsO.indexOf(gen2) < famsO.indexOf(gen3));
    assertTrue(MultiFamilyOrdering.isMonogamous(famsO));
    assertEquals(1, gen1.getFatherDistinctMates());
    assertEquals(1, gen1.getMotherDistinctMates());
    assertEquals(1, gen2.getFatherDistinctMates());
    assertEquals(1, gen2.getMotherDistinctMates());
    assertEquals(1, gen1b.getFatherDistinctMates());
    assertEquals(1, gen1b.getMotherDistinctMates());
    assertEquals(1, gen2b.getFatherDistinctMates());
    assertEquals(1, gen2b.getMotherDistinctMates());
    assertEquals(1, gen3.getFatherDistinctMates());
    assertEquals(1, gen3.getMotherDistinctMates());
  }

  public void testUnfaithful() {
    final Set<Family> fams = new LinkedHashSet<>();
    final Family gen1 = new Family("father", "mother", "child1", "child2", "child3");
    final Family gen1b = new Family("father", "newWife", "child4", "child5");
    fams.add(gen1b);
    fams.add(gen1);
    final List<Family> famsO = MultiFamilyOrdering.orderFamiliesAndSetMates(fams);
    assertEquals(2, famsO.size());
    assertTrue(famsO.containsAll(Arrays.asList(gen1, gen1b)));
    assertFalse(MultiFamilyOrdering.isMonogamous(famsO));
    assertEquals(2, gen1.getFatherDistinctMates());
    assertEquals(1, gen1.getMotherDistinctMates());
    assertEquals(2, gen1b.getFatherDistinctMates());
    assertEquals(1, gen1b.getMotherDistinctMates());
  }
  public void testUnfaithfulWithGKids() {
    final Set<Family> fams = new LinkedHashSet<>();
    final Family gen1 = new Family("father", "mother", "child1", "child2", "child3");
    final Family gen1b = new Family("father", "newWife", "child4", "child5");
    final Family gen2 = new Family("child1", "dInLaw1", "gchild1");
    final Family gen2b = new Family("sInLaw", "child4", "gchild2");
    fams.add(gen1b);
    fams.add(gen2);
    fams.add(gen2b);
    fams.add(gen1);
    final List<Family> famsO = MultiFamilyOrdering.orderFamiliesAndSetMates(fams);
    assertEquals(4, famsO.size());
    assertTrue(famsO.containsAll(Arrays.asList(gen1, gen1b, gen2, gen2b)));
    assertTrue(famsO.indexOf(gen1) < famsO.indexOf(gen2));
    assertTrue(famsO.indexOf(gen1) < famsO.indexOf(gen2b));
    assertTrue(famsO.indexOf(gen1b) < famsO.indexOf(gen2));
    assertTrue(famsO.indexOf(gen1b) < famsO.indexOf(gen2b));
    assertFalse(MultiFamilyOrdering.isMonogamous(famsO));
    assertEquals(2, gen1.getFatherDistinctMates());
    assertEquals(1, gen1.getMotherDistinctMates());
    assertEquals(2, gen1b.getFatherDistinctMates());
    assertEquals(1, gen1b.getMotherDistinctMates());
    assertEquals(1, gen2.getFatherDistinctMates());
    assertEquals(1, gen2.getMotherDistinctMates());
    assertEquals(1, gen2b.getFatherDistinctMates());
    assertEquals(1, gen2b.getMotherDistinctMates());
  }

  public void testUnfaithfulWithIncestuousHalfSibs() {
    final Set<Family> fams = new LinkedHashSet<>();
    final Family gen1 = new Family("father", "mother", "child1", "child2", "child3");
    final Family gen1b = new Family("father", "newWife", "child4", "child5");
    final Family gen2 = new Family("child1", "child4", "gchild1");
    fams.add(gen1b);
    fams.add(gen2);
    fams.add(gen1);
    final List<Family> famsO = MultiFamilyOrdering.orderFamiliesAndSetMates(fams);
    assertEquals(3, famsO.size());
    assertTrue(famsO.containsAll(Arrays.asList(gen1, gen1b, gen2)));
    assertTrue(famsO.indexOf(gen1) < famsO.indexOf(gen2));
    assertTrue(famsO.indexOf(gen1b) < famsO.indexOf(gen2));
    assertFalse(MultiFamilyOrdering.isMonogamous(famsO));

    assertEquals(2, gen1.getFatherDistinctMates());
    assertEquals(1, gen1.getMotherDistinctMates());
    assertEquals(2, gen1b.getFatherDistinctMates());
    assertEquals(1, gen1b.getMotherDistinctMates());
    assertEquals(1, gen2.getFatherDistinctMates());
    assertEquals(1, gen2.getMotherDistinctMates());
  }

  public void testBigDaddy() {
    final Set<Family> fams = new LinkedHashSet<>();
    final Family gen1 = new Family("father", "mother", "child1", "child2", "child3");
    final Family gen2 = new Family("father", "child2", "gchild1", "gchild2");
    final Family gen3 = new Family("father", "gchild1", "ggchild1", "ggchild2");
    final Family gen4 = new Family("father", "ggchild2", "gggchild1", "gggchild2");
    fams.add(gen2);
    fams.add(gen3);
    fams.add(gen4);
    fams.add(gen1);
    final List<Family> famsO = MultiFamilyOrdering.orderFamiliesAndSetMates(fams);
    assertEquals(4, famsO.size());
    assertTrue(famsO.containsAll(Arrays.asList(gen1, gen2, gen3, gen4)));
    assertTrue(famsO.indexOf(gen1) < famsO.indexOf(gen2));
    assertTrue(famsO.indexOf(gen2) < famsO.indexOf(gen3));
    assertTrue(famsO.indexOf(gen3) < famsO.indexOf(gen4));
    assertFalse(MultiFamilyOrdering.isMonogamous(famsO));
    assertEquals(4, gen1.getFatherDistinctMates());
    assertEquals(1, gen1.getMotherDistinctMates());
    assertEquals(4, gen2.getFatherDistinctMates());
    assertEquals(1, gen2.getMotherDistinctMates());
    assertEquals(4, gen3.getFatherDistinctMates());
    assertEquals(1, gen3.getMotherDistinctMates());
    assertEquals(4, gen4.getFatherDistinctMates());
    assertEquals(1, gen4.getMotherDistinctMates());
  }

  public void testBigDaddyNewWife() {
    final Set<Family> fams = new LinkedHashSet<>();
    final Family gen1 = new Family("father", "mother", "child1", "child2", "child3");
    final Family gen1b = new Family("father", "newWife", "child4", "child5");
    final Family gen2 = new Family("father", "child3", "gchild1");
    final Family gen2b = new Family("father", "child4", "gchild2");
    fams.add(gen2b);
    fams.add(gen1b);
    fams.add(gen2);
    fams.add(gen1);
    final List<Family> famsO = MultiFamilyOrdering.orderFamiliesAndSetMates(fams);
    assertEquals(4, famsO.size());
    assertTrue(famsO.containsAll(Arrays.asList(gen1, gen1b, gen2, gen2b)));
    assertTrue(famsO.indexOf(gen1) < famsO.indexOf(gen2));
    assertTrue(famsO.indexOf(gen1b) < famsO.indexOf(gen2b));
    assertFalse(MultiFamilyOrdering.isMonogamous(famsO));

    assertEquals(4, gen1.getFatherDistinctMates());
    assertEquals(1, gen1.getMotherDistinctMates());
    assertEquals(4, gen1b.getFatherDistinctMates());
    assertEquals(1, gen1b.getMotherDistinctMates());
    assertEquals(4, gen2.getFatherDistinctMates());
    assertEquals(1, gen2.getMotherDistinctMates());
    assertEquals(4, gen2b.getFatherDistinctMates());
    assertEquals(1, gen2b.getMotherDistinctMates());
  }

  public void testBigDaddyGchild() {
    final Set<Family> fams = new LinkedHashSet<>();
    final Family gen1 = new Family("father", "mother", "child1", "child2", "child3");
    final Family gen1b = new Family("father", "newWife", "child4", "child5");
    final Family gen2 = new Family("child1", "child4", "gchild1");
    final Family gen3 = new Family("father", "gchild1", "ggchild2");
    fams.add(gen3);
    fams.add(gen1b);
    fams.add(gen2);
    fams.add(gen1);
    final List<Family> famsO = MultiFamilyOrdering.orderFamiliesAndSetMates(fams);
    assertEquals(4, famsO.size());
    assertTrue(famsO.containsAll(Arrays.asList(gen1, gen1b, gen2, gen3)));
    assertTrue(famsO.indexOf(gen1) < famsO.indexOf(gen2));
    assertTrue(famsO.indexOf(gen1b) < famsO.indexOf(gen2));
    assertTrue(famsO.indexOf(gen2) < famsO.indexOf(gen3));
    assertFalse(MultiFamilyOrdering.isMonogamous(famsO));

    assertEquals(3, gen1.getFatherDistinctMates());
    assertEquals(1, gen1.getMotherDistinctMates());
    assertEquals(3, gen1b.getFatherDistinctMates());
    assertEquals(1, gen1b.getMotherDistinctMates());
    assertEquals(1, gen2.getFatherDistinctMates());
    assertEquals(1, gen2.getMotherDistinctMates());
    assertEquals(3, gen3.getFatherDistinctMates());
    assertEquals(1, gen3.getMotherDistinctMates());
  }

  public void testBigDaddyAndBigMommaGchild() {
    final Set<Family> fams = new LinkedHashSet<>();
    final Family gen1 = new Family("father", "mother", "child1", "child2", "child3");
    final Family gen1b = new Family("father", "newWife", "child4", "child5");
    final Family gen2 = new Family("child1", "mother", "gchild1");
    final Family gen3 = new Family("father", "gchild1", "ggchild2");
    fams.add(gen3);
    fams.add(gen1b);
    fams.add(gen2);
    fams.add(gen1);
    final List<Family> famsO = MultiFamilyOrdering.orderFamiliesAndSetMates(fams);
    assertEquals(4, famsO.size());
    assertTrue(famsO.containsAll(Arrays.asList(gen1, gen1b, gen2, gen3)));
    assertTrue(famsO.indexOf(gen1) < famsO.indexOf(gen2));
    assertTrue(famsO.indexOf(gen1b) < famsO.indexOf(gen2));
    assertTrue(famsO.indexOf(gen2) < famsO.indexOf(gen3));
    assertFalse(MultiFamilyOrdering.isMonogamous(famsO));

    assertEquals(3, gen1.getFatherDistinctMates());
    assertEquals(2, gen1.getMotherDistinctMates());
    assertEquals(3, gen1b.getFatherDistinctMates());
    assertEquals(1, gen1b.getMotherDistinctMates());
    assertEquals(1, gen2.getFatherDistinctMates());
    assertEquals(2, gen2.getMotherDistinctMates());
    assertEquals(3, gen3.getFatherDistinctMates());
    assertEquals(1, gen3.getMotherDistinctMates());
  }
}
