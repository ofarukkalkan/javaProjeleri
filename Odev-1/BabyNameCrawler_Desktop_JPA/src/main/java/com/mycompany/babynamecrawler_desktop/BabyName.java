/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.babynamecrawler_desktop;

import java.io.Serializable;
import javax.persistence.Column;
import javax.persistence.EmbeddedId;
import javax.persistence.Entity;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.Table;
import javax.xml.bind.annotation.XmlRootElement;

/**
 *
 * @author gldev
 */
@Entity
@Table(name = "BABYNAME")
@XmlRootElement
@NamedQueries({
    @NamedQuery(name = "BabyName.findAll", query = "SELECT b FROM BabyName b")
    , @NamedQuery(name = "BabyName.findByYearr", query = "SELECT b FROM BabyName b WHERE b.babyNamePK.yearr = :yearr")
    , @NamedQuery(name = "BabyName.findByName", query = "SELECT b FROM BabyName b WHERE b.babyNamePK.name = :name")
    , @NamedQuery(name = "BabyName.findByGender", query = "SELECT b FROM BabyName b WHERE b.babyNamePK.gender = :gender")
    , @NamedQuery(name = "BabyName.findByCount", query = "SELECT b FROM BabyName b WHERE b.count = :count")})
public class BabyName implements Serializable {

    private static final long serialVersionUID = 1L;
    @EmbeddedId
    protected BabyNamePK babyNamePK;
    @Column(name = "COUNT")
    private Integer count;

    public BabyName() {
    }

    public BabyName(BabyNamePK babyNamePK) {
        this.babyNamePK = babyNamePK;
    }

    public BabyName(int yearr, String name, Character gender) {
        this.babyNamePK = new BabyNamePK(yearr, name, gender);
    }

    public BabyNamePK getBabyNamePK() {
        return babyNamePK;
    }

    public void setBabyNamePK(BabyNamePK babyNamePK) {
        this.babyNamePK = babyNamePK;
    }

    public Integer getCount() {
        return count;
    }

    public void setCount(Integer count) {
        this.count = count;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (babyNamePK != null ? babyNamePK.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof BabyName)) {
            return false;
        }
        BabyName other = (BabyName) object;
        if ((this.babyNamePK == null && other.babyNamePK != null) || (this.babyNamePK != null && !this.babyNamePK.equals(other.babyNamePK))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "com.mycompany.babynamecrawler_desktop.BabyName[ babyNamePK=" + babyNamePK + " ]";
    }
    
}
