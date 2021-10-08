using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class AttackParticleController : MonoBehaviour
{
    private void FinishAttack()
    {
        Destroy(gameObject);
    }
}
