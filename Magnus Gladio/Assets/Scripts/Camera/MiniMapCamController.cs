using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MiniMapCamController : MonoBehaviour
{
    public float followSpeed = 2f;
    public float yOffSet = 2.5f;
    public float xOffSet = 0.5f;

  
    public Transform target;

    

    void Start()
    {
        SetupCamera();
    }

    public void SetupCamera()
    {
        target = GameObject.FindGameObjectWithTag("Player").GetComponent<Transform>();

    }
   
    // Update is called once per frame
    void LateUpdate()
    {
        if (target != null)
        {
            Vector3 newPosition = target.position;
            newPosition.x += xOffSet * target.GetComponent<PlayerController>().GetMovementDirectionOfPlayer();
            newPosition.y += yOffSet;
            newPosition.z = -10;
            transform.position = Vector3.Slerp(transform.position, newPosition, followSpeed * Time.deltaTime);
          
        }
    }

}
