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
