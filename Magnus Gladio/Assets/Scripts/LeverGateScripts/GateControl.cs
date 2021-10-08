using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GateControl : MonoBehaviour
{
    private GameObject cam;

    public void gateCam(GameObject gateCam)
    {
        cam = gateCam;
    }

    public void DisableCam()
    {
        cam.SetActive(false);
    }

}
    