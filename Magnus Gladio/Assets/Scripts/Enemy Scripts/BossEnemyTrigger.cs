using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class BossEnemyTrigger : MonoBehaviour
{

    public GameObject BossGate;
    public Transform bossCameraPoint;
    public Transform gateTargetPos;
    public float gateSpeed;

    private GameObject camObj;

    private bool 
        gateClosed = false, 
        camMoved = false;

    private void Awake()
    {
        camObj = GameObject.FindGameObjectWithTag("MainCamera");
    }

    public void OnTriggerEnter2D(Collider2D collision)
    {
        if (collision.CompareTag("Player"))
        {
            if (!gateClosed && !camMoved)
            {
                StartCoroutine(CloseGate());
                StartCoroutine(MoveCamera());
            }
        }
    }


    IEnumerator CloseGate()
    {
        gateClosed = true;
        yield return null;
        float t = 0;

        float d = Vector3.Distance(BossGate.transform.position, gateTargetPos.position);


        while (t <= 1 && Vector3.Distance(BossGate.transform.position, gateTargetPos.position) > 0.5f)
        {
            BossGate.transform.position = Vector3.Lerp(BossGate.transform.position, gateTargetPos.position,t);

            t += Time.deltaTime;
        }
        
    }

    IEnumerator MoveCamera()
    {
        camMoved = true;
        yield return null;
        float t = 0 ;
        while ( t <= 1 && Vector3.Distance(camObj.transform.position, bossCameraPoint.position) > 0.5f)
        {
            camObj.transform.position = Vector3.Lerp(camObj.transform.position, bossCameraPoint.position, t);
            t += Time.deltaTime;
        }
        

    }


}
