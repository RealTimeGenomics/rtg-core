/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.alignment;

import java.io.IOException;

import com.rtg.mode.Protein;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.protein.GotohProteinEditDistance;
import com.rtg.util.InvalidParamsException;

import junit.framework.TestCase;


/**
 */
public class ProteinActionsValidatorTest extends TestCase {

  public void testProteinEquals() throws InvalidParamsException, IOException {
    checkIsValidProtein(
        "ARNDCQEGHILKMFPSTWYV",
        "ARNDCQEGHILKMFPSTWYV"
        );
  }

  public void testProteinSubs() throws InvalidParamsException, IOException {
    checkIsValidProtein(
        "ARNDCQEGHILKMFPSTWYV",
        "ARNACQEGHLIKMGPSTWYV"
        );
  }

  public void testProteinDelete() throws InvalidParamsException, IOException {
    checkIsValidProtein(
        "ARNDCQEGHIIILLLLKMFPSTWYV",
        "ARNACQEGHI" + "LKMGPSTWYV"
        );
  }

  public void testProteinInserts() throws InvalidParamsException, IOException {
    checkIsValidProtein(
        "ARNDCQE    GHILKMF    PSTWYV".replaceAll(" ", ""),
         "RNACQEAVAVGHILKMFFFFFPSTWY"
        );
  }

  /** An example where ProteinEditDistance returns a non-optimal alignment. */
  public void testProteinBadAlignment() throws InvalidParamsException, IOException {
    checkIsValidProtein(
        "CEQGIHIILLLLKMFPSTWYV",
        "CQEGHI" + "LKMGPSTWYV"
        );
  }

  public void checkIsValidProtein(final String readString, final String templateString) throws InvalidParamsException, IOException {
    checkIsValid(readString, templateString, new ProteinScoringMatrix());
  }

  /**
   * Test the ActionsValidator by mutating entries in the actions array.
   * @param readString the read
   * @param templateString the template
   * @param matrix a ProteinScoringMatrix, or an Integer (for the DNA gap open penalty).
   */
  public void checkIsValid(final String readString, final String templateString, final ProteinScoringMatrix matrix) {
    final ActionsValidator validator = new ProteinActionsValidator(matrix);
    final byte[] read = Protein.encodeProteins(readString);
    final byte[] temp = Protein.encodeProteins(templateString);
    final int maxScore = Integer.MAX_VALUE;
    final int msf = Math.max((int) (read.length * 0.7), 7);

    int[] a0 = new GotohProteinEditDistance(matrix).calculateEditDistance(read, read.length, temp, 0, maxScore, msf, true);
    a0 = ActionsHelper.copy(a0); // make it minimal.
    assertTrue(validator.isValid(a0, read, read.length, temp, false, maxScore));

    // Now try perturbing each integer in the actions array,
    // and check that every change gives a false result from isValid.
    final int len = ActionsHelper.actionsCount(a0);
    for (int i = 0; i <= ActionsHelper.ACTIONS_START_INDEX + (len - 1) / ActionsHelper.ACTIONS_PER_INT; ++i) {
      //System.err.println("\nPerturbing entry a0[" + i + "] = " + a0[i]);
      a0[i]--;
      assertFalse(validator.mErrorMsg, validator.isValid(a0, read, read.length, temp, false, maxScore));
      a0[i]++;
      assertTrue(validator.mErrorMsg, validator.isValid(a0, read, read.length, temp, false, maxScore));
      a0[i]++;
      assertFalse(validator.isValid(a0, read, read.length, temp, false, maxScore));
      a0[i]--;
      assertTrue(validator.isValid(a0, read, read.length, temp, false, maxScore));
    }
  }
}
