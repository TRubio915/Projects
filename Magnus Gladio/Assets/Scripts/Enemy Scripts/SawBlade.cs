using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SawBlade : MonoBehaviour
{
    public float rotationSpeed; //Speed of blade

    void Start()
    {
        rotationSpeed = 5;
    }

    void FixedUpdate()
    {
        this.transform.Rotate(new Vector3 (0, 0, rotationSpeed));
    }
}
