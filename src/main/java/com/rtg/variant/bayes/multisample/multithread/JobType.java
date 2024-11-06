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

package com.rtg.variant.bayes.multisample.multithread;

import java.util.List;

import com.rtg.scheduler.Result;
import com.rtg.variant.bayes.multisample.Complexities;

/**
 * Execution phases within one chunk.
 */
public enum JobType {
  /** Increment model counts and do initial calling. */
  INCR {
    @Override
    public boolean validArguments(Result[] args) {
      return args.length == 0;
    }
    @Override
    public boolean validResult(Result result) {
      return result.length() == 2
          && (result.result(0) instanceof Complexities)
          && (result.result(1) instanceof Integer)
          ;
    }
  },

  /** Tidy dangling regions at the ends of chunks. */
  DANGLING {
    @Override
    public boolean validArguments(Result[] args) {
      if (args.length != 3) {
        return false;
      }
      if (args[0] == null) {
        return INCR.validResult(args[1]);
      }
      if (args[1] == null) {
        return INCR.validResult(args[0]);
      }
      return INCR.validResult(args[0]) && INCR.validResult(args[1]);
    }
    @Override
    public boolean validResult(Result result) {
      return result.length() == 1
          && (result.result(0) == null || result.result(0) instanceof Complexities)
          ;
    }
  },

  /** Complex calling. */
  COMPLEX {
    @Override
    public boolean validArguments(Result[] args) {
      if (args.length != 1) {
        return false;
      }
      return DANGLING.validResult(args[0]);
    }
    @Override
    public boolean validResult(Result result) {
      return result.length() == 2
        && (result.result(0) == null || result.result(0) instanceof List)
        && (result.result(1) == null || result.result(1) instanceof Complexities)
        ;
    }
  },
  /** flushes the circular buffer */
  FLUSH {
    @Override
    public boolean validArguments(Result[] args) {
      if (args.length != 2) {
        return false;
      }
      if (args[0] == null) {
        return COMPLEX.validResult(args[1]);
      }
      if (args[1] == null) {
        return INCR.validResult(args[0]);
      }
      return INCR.validResult(args[0]) && COMPLEX.validResult(args[1]);
    }

    @Override
    public boolean validResult(Result result) {
      return result.length() == 0;
    }
  },
  /** Merge results from initial calling and complex calling. */
  FILTER {
    @Override
    public boolean validArguments(Result[] args) {
      if (args.length != 2) {
        return false;
      }
      return (args[0] == null || INCR.validResult(args[0])) && (args[1] == null || COMPLEX.validResult(args[1]));
    }
    @Override
    public boolean validResult(Result result) {
      return result.length() == 2
          && (result.result(0) == null || result.result(0) instanceof List)
          && (result.result(1) == null || result.result(1) instanceof Integer)
          ;
    }
  },
  /** Write regions to <code>regions.bed</code> file. */
  BED {
    @Override
    public boolean validArguments(Result[] args) {
      if (args.length != 2) {
        return false;
      }
      if (args[0] == null) {
        //can happen at time 0
        return true;
      }
      return BED.validResult(args[0]) && COMPLEX.validResult(args[1]);
    }
    @Override
    public boolean validResult(Result result) {
      return result.length() == 0;
    }
  },
  /** Write called Variants. */
  OUT {
    @Override
    public boolean validArguments(Result[] args) {
      return args.length == 3
        && (args[0] == null || OUT.validResult(args[0]))
        && (args[1] == null || FILTER.validResult(args[1]))
        && (args[2] == null || FLUSH.validResult(args[2]));
    }
    @Override
    public boolean validResult(Result result) {
      return result == null || (result.length() == 1 && (result.result(0) == null || result.result(0) instanceof List));
    }
  };

  /**
   * @param args to be validated.
   * @return true iff the number and types of the results are correct as arguments for a job with this identifier.
   */
  public abstract boolean validArguments(Result[] args);

  /**
   * @param result to be validated.
   * @return true iff the number and types of the objects in result are correct as the result for a job with this identifier.
   */
  public abstract boolean validResult(Result result);


}
