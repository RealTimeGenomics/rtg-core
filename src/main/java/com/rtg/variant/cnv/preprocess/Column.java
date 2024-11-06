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
package com.rtg.variant.cnv.preprocess;

import com.reeltwo.jumble.annotations.TestClass;

/**
 * Holds a column of data
 */
@TestClass("com.rtg.variant.cnv.preprocess.StringColumnTest")
public abstract class Column implements Cloneable {

  private String mName;

  Column(String name) {
    mName = name;
  }

  public void setName(String name) {
    mName = name;
  }

  /**
   * @return the name of the column
   */
  public String getName() {
    return mName;
  }

  /**
   * @return the number of items
   */
  public abstract int size();

  /**
   * Add a value
   * @param value the value to add
   */
  abstract void add(String value);

  /**
   * Add a value at a specific index
   * @param i the index of the entry to insert at
   * @param value the value to add
   */
  abstract void add(int i, String value);

  /**
   * Removes the element at the specified position in this list.
   * Shifts any subsequent elements to the left (subtracts one from their indices).
   * @param i the index of the entry to remove
   */
  abstract void remove(int i);

  /**
   * Get a value in String representation
   * @param i the item
   * @return the String representation of the item
   */
  abstract String toString(int i);

  abstract Column filter(RegionPredicate p);
}
