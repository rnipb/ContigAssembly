package com.rnipb.pab.alignment;

import java.io.Serializable;

public class Pair<T1, T2> implements Serializable {
    private T1 first;
    private T2 second;

    public T1 first() {
        return this.first;
    }

    public T2 second() {
        return this.second;
    }

    public Pair(T1 first, T2 second) {
        this.first = first;
        this.second = second;
    }

    public String toString() {
        return "(" + this.first + "," + this.second + ")";
    }

    public boolean equals(Object o) {
        if (this == o) {
            return true;
        } else if (o != null && this.getClass() == o.getClass()) {
            Pair pair = (Pair)o;
            if (this.first != null) {
                if (!this.first.equals(pair.first)) {
                    return false;
                }
            } else if (pair.first != null) {
                return false;
            }

            if (this.second != null) {
                if (this.second.equals(pair.second)) {
                    return true;
                }
            } else if (pair.second == null) {
                return true;
            }

            return false;
        } else {
            return false;
        }
    }

    public int hashCode() {
        int result = this.first != null ? this.first.hashCode() : 0;
        result = 31 * result + (this.second != null ? this.second.hashCode() : 0);
        return result;
    }
}