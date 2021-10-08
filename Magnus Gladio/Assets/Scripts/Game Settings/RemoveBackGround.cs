using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RemoveBackGround : MonoBehaviour
{

    public GameObject backGround;

    private void OnTriggerEnter2D(Collider2D collision)
    {
        if(collision.gameObject.CompareTag("Player"))
        {
            backGround.SetActive(false);
        }
    }

   
}
