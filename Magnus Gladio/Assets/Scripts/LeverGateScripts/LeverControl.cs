using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class LeverControl : MonoBehaviour
{


    public GameObject gateLinked;

    private Animator anim;
    private Collider2D coll;

    public GameObject gateCamUI;

    public GameObject gateCam;

    public bool enableLever = false;

    // Start is called before the first frame update
    void Start()
    {
        anim = GetComponent<Animator>();
        coll = GetComponent<BoxCollider2D>();
    }
    private void Update()
    {
        if(enableLever)
        {
            TurnLever();
            enableLever = false;
        }

        if(Input.GetKeyDown(KeyCode.Alpha4))
        {
            TurnLever();
        }
    }

    public void TurnLever()
    {
        anim.enabled = true;
        coll.enabled = false;
        gateCam.transform.position = new Vector3(gateLinked.transform.localPosition.x, gateLinked.transform.localPosition.y,gateCam.transform.position.z);
        gateCamUI.SetActive(true);
    }

    public void LeverOn()
    {
        gateLinked.GetComponent<Animator>().enabled = true;
        gateLinked.SendMessage("gateCam", gateCamUI);
    }

}
