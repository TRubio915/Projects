using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class DropAbility : MonoBehaviour
{


    public GameObject abilityPrefab;

    
    public void AbilityDropOnDeath()
    {
        if (abilityPrefab != null)
        {
            GameObject abilityClone = Instantiate(abilityPrefab, transform.position, Quaternion.identity);

            abilityClone.GetComponent<Rigidbody2D>().AddForce(Vector2.up * 50);
        }
    }


}
