using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class BossFinal : Enemy
{
    public Transform target; //Target of enemy, typically the player
    public float aggroRange; //Aggro Range for enemy
    void Start()
    {
        target = GameObject.FindGameObjectWithTag("Player").transform;
    }

    void Update()
    {
        
    }

    private void OnDrawGizmos()
    {
        base.OnDrawGizmos();
        //Draw lines so detection is more clear
        Gizmos.DrawWireSphere(transform.position, aggroRange);
        Gizmos.color = Color.white;
    }
}
