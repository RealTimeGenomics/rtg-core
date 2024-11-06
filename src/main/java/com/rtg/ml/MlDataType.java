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
package com.rtg.ml;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

/**
 * Handle saving and loading of machine learning data values, may be extended later to provide other functions
 * (converting, splitting...)
 */
public enum MlDataType {

  /** values are strings */
  STRING {
    @Override
    public void save(Object value, DataOutputStream dos) throws IOException {
      dos.writeUTF((String) value);
    }

    @Override
    public Object load(DataInputStream dis) throws IOException {
      return dis.readUTF();
    }

    @Override
    public boolean isNumeric() {
      return false;
    }

  },
  /** values are boolean */
  BOOLEAN {
    @Override
    public void save(Object value, DataOutputStream dos) throws IOException {
      dos.writeBoolean((Boolean) value);
    }

    @Override
    public Object load(DataInputStream dis) throws IOException {
      return dis.readBoolean();
    }

    @Override
    public boolean isNumeric() {
      return false;
    }

  },
  /** values are integer */
  INTEGER {
    @Override
    public void save(Object value, DataOutputStream dos) throws IOException {
      dos.writeInt((Integer) value);
    }

    @Override
    public Object load(DataInputStream dis) throws IOException {
      return dis.readInt();
    }

    @Override
    public boolean isNumeric() {
      return true;
    }

  },
  /** values are double */
  DOUBLE {
    @Override
    public void save(Object value, DataOutputStream dos) throws IOException {
      dos.writeDouble((Double) value);
    }

    @Override
    public Object load(DataInputStream dis) throws IOException {
      return dis.readDouble();
    }

    @Override
    public boolean isNumeric() {
      return true;
    }

  };

  /**
   * saves an object of the corresponding data type to the stream
   * @param value the value to save
   * @param dos the stream to save to
   * @throws IOException if an IO error occurs
   */
  public abstract void save(Object value, DataOutputStream dos) throws IOException;

  /**
   * loads an object of the corresponding data type from the stream
   * @param dis the stream to load from
   * @return the value
   * @throws IOException if an IO error occurs
   */
  public abstract Object load(DataInputStream dis) throws IOException;

  /***
   * @return true if type is numeric (i.e. values are ordered and can be interpolated)
   */
  public abstract boolean isNumeric();

}
