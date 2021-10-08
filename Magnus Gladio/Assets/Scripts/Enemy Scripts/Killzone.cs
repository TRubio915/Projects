using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Killzone : MonoBehaviour
{
    public float damage;
    public void OnCollisionEnter2D(Collision2D other)
    {
        if(other.gameObject.tag.Equals("Player"))
        {
            float[] attackDetails = new float[2];
            
            attackDetails[0] = damage;
            
            other.transform.SendMessage("Damage", attackDetails);
            
        }
    }
}
