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
 * Generate a single homozygous mutation on a template using a specification. The specification consists of
 * a string with the following grammar:
 * <pre>
 * single ::= {[count]{"X", "I", "J", "D", "=", "E"}}+
 *
 * The codes specify:
 * X mutate to a valid nucleotide not equal to the template but may be equal to a previous mutation.
 * Y mutate to a valid nucleotide not equal to the template nor equal to a previous mutation.
 * I insert a random valid nucleotide into the mutated result.
 * J insert a random valid nucleotide into the mutated result not equal to a previous mutation at that point in the mutated result.
 * D delete a nucleotide from the template.
 * = or E use the same nucleotide as the template.
 * </pre>
 * The class is re-entrant and thread safe.
 */
class MutatorSingle extends IntegralAbstract {

  private final String mSpecification;

  private final int mRefLength;
  private final int mMutLength;

  /**
   * Expand a spec with counts in it e.g. <code>5I</code> becomes <code>IIIII</code>
   * @param spec containing counts e.g. <code>5I</code>
   * @return a spec with the counted mutations expanded
   */
  static String expandSpec(String spec) {
    final StringBuilder sb = new StringBuilder();
    int count = 0;
    for (int s = 0; s < spec.length(); ++s) {
      final char c = spec.charAt(s);
      if (c >= '0' && c <= '9') {
        count = count * 10 + c - '0';
      } else {
        if (count == 0) {
          count = 1;
        }
        for (int i = 0; i < count; ++i) {
          sb.append(c);
        }
        count = 0;
      }
    }
    if (count != 0) {
      throw new RuntimeException(spec);
    }
    return sb.toString();
  }

  /**
   * @param specification for the mutation. See class documentation.
   */
  MutatorSingle(String specification) {
    mSpecification = expandSpec(specification);
    int t = 0;
    int r = 0;
    for (int s = 0; s < mSpecification.length(); ++s) {
      final char c = mSpecification.charAt(s);
      switch (c) {
        case 'X':
        case 'Y':
        case '=':
        case 'E':
          ++t;
          ++r;
          break;
        case 'I':
        case 'J':
          ++r;
          break;
        case 'D':
          ++t;
          break;
        default:
          throw new RuntimeException(mSpecification);
      }
    }
    mRefLength = t;
    mMutLength = r;
  }

  /**
   * @return length of mutation on reference
   */
  public int getReferenceLength() {
    return mRefLength;
  }

  /**
   * @return length of mutated result
   */
  public int getMutationLength() {
    return mMutLength;
  }

  /**
   * Almost but not quite the same as the method in the Mutator interface. The difference is that it checks another
   * result to make sure the mutation is different (for the case when this is the second of a heterozygous pair).
   * @param template the template to mutate
   * @param position the zero based position in the template to generate a mutation
   * @param random number generator used for generating mutations.
   * @param other a previous mutation - a check is made that an X code differs from this - will be ignored if it is null.
   * @return the generated mutation
   */
  MutatorResult generateMutation(byte[] template, int position, PortableRandom random, byte[] other) {
    //guaranteed to be long enough to hold the result
    //doing this here rather than as a field ensures the class is re-entrant.
    final byte[] result = new byte[mSpecification.length()];
    int t = position;
    int r = 0;
    for (int s = 0; s < mSpecification.length(); ++s) {
      assert 0 <= r && r <= result.length;
      assert position <= t && t <= template.length;
      final char c = mSpecification.charAt(s);
      switch (c) {
        case 'X':
          if (t >= template.length) {
            return null;
          }
          result[r] = minus(random, template[t]);
          ++t;
          ++r;
          break;
        case 'Y':
          if (t >= template.length) {
            return null;
          }
          final byte tem = template[t];
          if (other == null) {
            result[r] = minus(random, tem);
          } else {
            result[r] = minus(random, tem, other[r]);
          }
          ++t;
          ++r;
          break;
        case 'I':
          result[r] = random(random);
          ++r;
          break;
        case 'J':
          if (other == null) {
            result[r] = random(random);
          } else {
            result[r] = minus(random, other[r]);
          }
          ++r;
          break;
        case 'D':
          ++t;
          break;
        case '=':
        case 'E':
          if (t >= template.length) {
            return null;
          }
          result[r] = template[t];
          ++t;
          ++r;
          break;
        default:
          throw new RuntimeException(mSpecification);
      }
    }
    final byte[] res = new byte[r];
    System.arraycopy(result, 0, res, 0, r);
    return new MutatorResult(res, res, t - position);
  }

  static byte minus(PortableRandom random, byte a, byte b) {
    if (a == b) {
      return minus(random, a);
    }
    if (a > b) {
      return minus(random, b, a);
    }
    final int ra = random.nextInt(2) + 1;
    final int sa = ra + (ra >= a ? 1 : 0);
    final int ta = sa + (sa >= b ? 1 : 0);
    return (byte) ta;
  }

  static byte minus(PortableRandom random, byte a) {
    final int ra = random.nextInt(3) + 1;
    final int sa = ra + (ra >= a ? 1 : 0);
    return (byte) sa;
  }

  static byte random(PortableRandom random) {
    return (byte) (random.nextInt(4) + 1);
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mSpecification);
    return true;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(mSpecification);
  }

}
