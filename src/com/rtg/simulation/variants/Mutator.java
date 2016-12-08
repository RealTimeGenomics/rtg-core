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

package com.rtg.simulation.variants;

import com.rtg.util.PortableRandom;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Mutate a template using a specification. The specification consists of
 * a string with the following grammar:
 * <pre>
 * single ::= {"X", "Y", "J", "I", "D", "=", "E"}+
 * specification ::= single [{":", "_" } single]
 *
 * The two halves of the specification apply to the two haplotypes (if there is only one then the
 * mutation is homozygous). The codes specify:
 * X mutate to a valid nucleotide not equal to the template but may be equal to a previous mutation.
 * Y mutate to a valid nucleotide not equal to the template nor equal to a previous mutation.
 * I insert a random valid nucleotide into the mutated result.
 * J insert a random valid nucleotide into the mutated result not equal to a previous mutation at that point in the mutated result.
 * D delete a nucleotide from the template.
 * = or E use the same nucleotide as the template.
 * </pre>
 */
public class Mutator extends IntegralAbstract {

  private final MutatorSingle mFirst;

  private final MutatorSingle mSecond;

  /**
   * @param specification single class docs
   */
  public Mutator(final String specification) {
    final String[] split = specification.split("_|:");
    switch (split.length) {
      case 1:
        mFirst = new MutatorSingle(split[0]);
        mSecond = null;
        break;
      case 2:
        int i = 0;
        for (; i < split[0].length() && i < split[1].length(); ++i) {
          if (split[0].charAt(split[0].length() - 1 - i) == 'E' && split[1].charAt(split[1].length() - 1 - i) == 'E') {
            continue;
          }
          break;
        }
        mFirst = new MutatorSingle(split[0].substring(0, split[0].length() - i));
        mSecond = new MutatorSingle(split[1].substring(0, split[1].length() - i));
        break;
      default:
        throw new RuntimeException(specification);
    }
    assert integrity();
  }

  /**
   * @return length of reference that is mutated
   */
  public int getReferenceLength() {
    assert mSecond == null || mSecond.getReferenceLength() == mFirst.getReferenceLength();
    return mFirst.getReferenceLength();
  }

  /**
   * @return true if mutation contains an insertion or deletion
   */
  public boolean isIndel() {
    return mFirst.getReferenceLength() != mFirst.getMutationLength() || (mSecond != null && mSecond.getReferenceLength() != mSecond.getMutationLength());
  }

  /**
   * Generates a mutation for the specified template.
   * @param template the template to mutate
   * @param position the zero based position in the template to generate a mutation
   * @param random number generator used for generating mutations.
   * @return the generated mutation (may be null if a valid mutation cannot be generated at that position).
   */
  public MutatorResult generateMutation(byte[] template, int position, PortableRandom random) {
    final MutatorResult first = mFirst.generateMutation(template, position, random, null);
    final byte[] firstHaplotype = first.getFirstHaplotype();
    if (mSecond == null) {
      return new MutatorResult(firstHaplotype, firstHaplotype, first.getConsumed());
    }
    final MutatorResult second = mSecond.generateMutation(template, position, random, firstHaplotype);
    assert first.getConsumed() == second.getConsumed();
    return new MutatorResult(firstHaplotype, second.getFirstHaplotype(), first.getConsumed());
  }

  @Override
  public void toString(StringBuilder sb) {
    final String right = mSecond != null ? ":" + mSecond : "";
    sb.append("Mutator:").append(mFirst).append(right);
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mFirst);
    return true;
  }

}
