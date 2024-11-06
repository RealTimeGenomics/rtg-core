/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.mode;

import static com.rtg.mode.SequenceMode.BIDIRECTIONAL;
import static com.rtg.mode.SequenceMode.PROTEIN;
import static com.rtg.mode.SequenceMode.TRANSLATED;
import static com.rtg.mode.SequenceMode.UNIDIRECTIONAL;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

import com.rtg.util.PseudoEnum;

/**
 * The translation modes possible with the various program modes.
 * Only this subset make any sense.
 */
public abstract class ProgramMode implements Serializable, PseudoEnum {

  /** slimn mode. */
  public static final ProgramMode
  SLIMN = new ProgramMode(UNIDIRECTIONAL, BIDIRECTIONAL, SequenceType.DNA, "SLIMN", 0) {
    @Override
    public ProgramMode flip() {
      return SLIMN_FLIP;
    }
  };

  /** slimp mode. */
  public static final ProgramMode
  SLIMP = new ProgramMode(PROTEIN, PROTEIN, SequenceType.PROTEIN, "SLIMP", 1) {
    @Override
    public ProgramMode flip() {
      return this;
    }
  };

  /** slimx mode. */
  public static final ProgramMode
  SLIMX = new ProgramMode(PROTEIN, TRANSLATED, SequenceType.PROTEIN, "SLIMX", 2) {
    @Override
    public ProgramMode flip() {
      return TSLIMN;
    }
  };

  /** tslimn mode. */
  public static final ProgramMode
  TSLIMN = new ProgramMode(TRANSLATED, PROTEIN, SequenceType.PROTEIN, "TSLIMN", 3) {
    @Override
    public ProgramMode flip() {
      return SLIMX;
    }
  };

  /** tslimx mode. */
  public static final ProgramMode
  TSLIMX = new ProgramMode(TRANSLATED, TRANSLATED, SequenceType.PROTEIN, "TSLIMX", 4) {
    @Override
    public ProgramMode flip() {
      return this;
    }
  };

  /** duster mode. */
  public static final ProgramMode
  DUSTER = new ProgramMode(UNIDIRECTIONAL, UNIDIRECTIONAL, SequenceType.DNA, "DUSTER", 5) {
    @Override
    public ProgramMode flip() {
      return this;
    }
  };

  /** Flip of slimn mode - only used for build on query never direct from user. */
  public static final ProgramMode
  SLIMN_FLIP = new ProgramMode(BIDIRECTIONAL, UNIDIRECTIONAL, SequenceType.DNA, "SLIMN_FLIP", 6) {
    @Override
    public ProgramMode flip() {
      return SLIMN;
    }
  };

  /** Only used for phylogeny computations - never direct from user. */
  public static final ProgramMode
  PHYLOGENY = new ProgramMode(UNIDIRECTIONAL, UNIDIRECTIONAL, SequenceType.DNA, "PHYLOGENY", 7) {
    @Override
    public ProgramMode flip() {
      return this;
    }
  };


  private static final Map<String, Object> VALUE_OF = new HashMap<>();

  private static final Object[] VALUES = {SLIMN, SLIMP, SLIMX, TSLIMN, TSLIMX};

  /**
   * Generate array of all the possible SequenceMode singletons.
   * These are in the same ordering as ordinal().
   * @return array of all the possible SequenceMode singletons.
   */
  public static Object[] values() {
    return VALUES.clone();
  }

  static {
    for (final Object sm : values()) {
      VALUE_OF.put(sm.toString(), sm);
    }
  }

  /**
   * Get the SequenceMode singleton with the specified value (aka name).
   * @param str the name of a SequenceMode singleton.
   * @return the singleton SequenceMode
   * @throws IllegalArgumentException if str is not a valid name.
   */
  public static ProgramMode valueOf(final String str) {
    final ProgramMode res = (ProgramMode) VALUE_OF.get(str);
    if (res == null) {
      throw new IllegalArgumentException(str);
    }
    return res;
  }

  private final SequenceMode mSubjectMode;

  private final SequenceMode mQueryMode;

  private final SequenceType mCodeType;

  private final String mToString;

  private final int mOrdinal;

  private ProgramMode(final SequenceMode subjectMode, final SequenceMode queryMode, final SequenceType codeType, final String toString, final int ordinal) {
    mSubjectMode = subjectMode;
    mQueryMode = queryMode;
    mCodeType = codeType;
    mToString = toString;
    mOrdinal = ordinal;
  }

  @Override
  public String name() {
    return toString();
  }

  /**
   * The sequence mode of the query sequences.
   * @return the sequence mode.
   */
  public SequenceMode queryMode() {
    return mQueryMode;
  }

  /**
   * The sequence mode of the subject sequences.
   * @return the sequence mode.
   */
  public SequenceMode subjectMode() {
    return mSubjectMode;
  }

  /**
   * Get the type of residue after translation (it is the same for query and subject).
   * @return the code type.
   */
  public SequenceType translatedType() {
    return mCodeType;
  }

  /**
   * Get the mode when subject and query are flipped.
   * @return the mode when subject and query are flipped.
   */
  public abstract ProgramMode flip();

  @Override
  public String toString() {
    return mToString;
  }

  /**
   * Get the unique integer ordinal value associated with each singleton.
   * @return the unique integer ordinal value associated with each singleton.
   */
  @Override
  public int ordinal() {
    return mOrdinal;
  }

  /**
   * Special handling of Serialization to ensure we get singletons on deserialization.
   * @return a singleton.
   */
  Object readResolve() {
    if (ordinal() == 5) {
      return DUSTER;
    }
    if (ordinal() == 6) {
      return SLIMN_FLIP;
    }
    return VALUES[ordinal()];
  }


}

