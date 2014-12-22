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
package com.rtg.mode;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.Locale;

import com.rtg.util.StringUtils;

import junit.framework.TestCase;


/**
 */
public class ProteinTest extends TestCase {

  /**
   * Test method for {@link com.rtg.mode.Protein()}.
   * @throws InvocationTargetException
   * @throws IllegalAccessException
   * @throws IllegalArgumentException
   * @throws NoSuchMethodException
   * @throws SecurityException
   */
  public final void test() throws IllegalArgumentException, IllegalAccessException, InvocationTargetException, SecurityException, NoSuchMethodException {
    // Check toString of values
    final String methodName = "values";
    Method m = Protein.class.getMethod(methodName);
    Protein[] r = (Protein[]) m.invoke(null);
    assertEquals("[X, *, A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V]", Arrays.toString(r));
    final String valueOfMethod = "valueOf";
    m = Protein.class.getMethod(valueOfMethod, String.class);
    // Check ordinal and valueOf
    for (int i = 0; i < r.length; i++) {
      assertEquals(i, r[i].ordinal());
      if (r[i] != Protein.STOP) {
        assertEquals(r[i], m.invoke(null, r[i].toString()));
        assertEquals(Protein.valueOf(r[i].toString()), r[i]);
      }
    }
  }

  //Copied  from the IUPAC page
  private static final String STR = ""
    + "X\tXaa\tAny amino acid" + StringUtils.LS
    + "*\t***\tTranslated stop codon" + StringUtils.LS
    + "A\tAla\tAlanine" + StringUtils.LS
    + "R\tArg\tArginine" + StringUtils.LS
    + "N\tAsn\tAsparagine" + StringUtils.LS
    + "D\tAsp\tAspartic acid" + StringUtils.LS
    + "C\tCys\tCysteine" + StringUtils.LS
    + "Q\tGln\tGlutamine" + StringUtils.LS
    + "E\tGlu\tGlutamic acid" + StringUtils.LS
    + "G\tGly\tGlycine" + StringUtils.LS
    + "H\tHis\tHistidine" + StringUtils.LS
    + "I\tIle\tIsoleucine" + StringUtils.LS
    + "L\tLeu\tLeucine" + StringUtils.LS
    + "K\tLys\tLysine" + StringUtils.LS
    + "M\tMet\tMethionine" + StringUtils.LS
    + "F\tPhe\tPhenylalanine" + StringUtils.LS
    + "P\tPro\tProline" + StringUtils.LS
    + "S\tSer\tSerine" + StringUtils.LS
    + "T\tThr\tThreonine" + StringUtils.LS
    + "W\tTrp\tTryptophan" + StringUtils.LS
    + "Y\tTyr\tTyrosine" + StringUtils.LS
    + "V\tVal\tValine" + StringUtils.LS;

  /**
   * Test method for {@link com.rtg.mode.Protein()}.
   */
  public final void testNames() {
    final StringBuilder sb = new StringBuilder();
    for (final Protein pr : Protein.values()) {
      sb.append(pr.name()).append("\t").append(pr.threeLetter()).append("\t").append(pr.fullName()).append(StringUtils.LS);
      assertEquals(0, pr.compareTo(pr));
      assertFalse(pr.equals(null));
      assertTrue(pr.equals(pr));
      assertEquals(pr.ordinal(), pr.hashCode());
    }
    assertFalse(Protein.A.equals(Protein.C));
    assertEquals(STR, sb.toString());
  }

  /**
   * Test method for {@link com.rtg.mode.Protein#ignore()}.
   */
  public final void testIgnore() {
    final Protein unknown = Protein.X;
    final Protein stop = Protein.STOP;
    assertTrue(unknown.ignore());
    assertEquals(0, unknown.ordinal());
    assertTrue(Protein.valueOf("X").ignore());
    for (final Protein protein : Protein.values()) {
      if (protein != unknown && protein != stop) {
        assertFalse(protein.ignore());
      }
      if (protein.ignore()) {
        assertTrue(protein.ordinal() < SequenceType.PROTEIN.firstValid());
      } else {
        assertTrue(protein.ordinal() >= SequenceType.PROTEIN.firstValid());
      }
    }
  }

  /**
   * Test method for {@link com.rtg.mode.Protein#type()}.
   */
  public final void testType() {
    for (final Protein protein : Protein.values()) {
      assertEquals(SequenceType.PROTEIN, protein.type());
    }
  }

  /**
   * Test method for {@link com.rtg.mode.Protein#toString()}.
   */
  public final void testToString() {
    for (final Protein protein : Protein.values()) {
      final String str = protein.toString();
      assertEquals(1, str.length());
      assertEquals(str.toUpperCase(Locale.getDefault()), str);
    }
  }

}

